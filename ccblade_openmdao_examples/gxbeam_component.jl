using ConcreteStructs
using ComponentArrays
using ForwardDiff
using GXBeam, LinearAlgebra
using OpenMDAO: AbstractExplicitComp

@concrete struct SolverComp <: AbstractExplicitComp
    rho
    E
    nu
    points
    apply_nonlinear_forwarddiffable!
    x
    y
    J
    forwarddiff_config
end

function create_assembly(; points, Tp, Np, omega, chord, twist, A, Iyy, Izz, Iyz, rho, E, nu)

    x = [points[i, 1][1] for i = 1:length(points)]
    nelems = length(points)-1
    start = 1:nelems
    stop = 2:nelems+1

    G = E/(2*(1 + nu))

    x3c = 0.0 *cos.(twist).*chord  # 0.17*cos.(twist).*chord
    ky = 10*(1 + nu)/(12 + 11*nu)
    kz = 10*(1 + nu)/(12 + 11*nu)

    H22 = E*Iyy
    H33 = E*Izz
    H23 = E*Iyz

    Ay = A./ky
    Az = A./kz

    stiffness = fill(zeros(6,6), nelems)
    mass = fill(zeros(6, 6), nelems)

    for i in 1:nelems
        elem_stiff = [[E*A[i], 0,       0,       0,  0,       0     ]#=
                    =#[0,      G*Ay[i], 0,       0,  0,       0     ]#=
                    =#[0,      0,       G*Az[i], 0,  0,       0     ]#=
                    =#[0,      0,       0,       0,  0,       0     ]#=
                    =#[0,      0,       0,       0,  H22[i], -H23[i]]#=
                    =#[0,      0,       0,       0, -H23[i],  H33[i]]]
        elem_mass = [[rho*A[i],        0,               0,        0,               rho*A[i]*x3c[i], 0     ]#=
                   =#[0,               rho*A[i],        0,       -rho*A[i]*x3c[i], 0,               0     ]#=
                   =#[0,               0,               rho*A[i], 0,               0,               0     ]#=
                   =#[0,              -rho*A[i]*x3c[i], 0,        Iyy[i]+Izz[i],   0,               0     ]#=
                   =#[rho*A[i]*x3c[i], 0,               0,        0,               Iyy[i],         -Iyz[i]]#=
                   =#[0,               0,               0,        0,              -Iyz[i],          Izz[i]]]
        stiffness[i] = elem_stiff
        mass[i] = elem_mass
    end

    # create assembly
    assembly = Assembly(points, start, stop, stiffness=stiffness, mass=mass)

    # set prescribed conditions (fixed left endpoint)
    bcs = Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
    )

    function interp_aero_force(x, radii, f_aero)
        # Linearly interpolate the aerodynamic loads
        # x: coordinate where we want to evaluate the force
        # radii: vector with coordinates corresponding to aero force vector
        # f_aero: vector with aero forces

        if x < radii[1] || x >= radii[end]  # x not in range of radii
            q = 0
        elseif x == radii[1]  # to prevent dividing by zero
            q = f_aero[1]
        else
            k = findfirst(radii .>= x)
            j = k - 1
            q = f_aero[j] + (x - radii[j])*(f_aero[k] - f_aero[j])/(radii[k] - radii[j])
        end

        return q
    end

    # define distributed loads
    distributed_loads = Dict()
    for i = 1:nelems
        distributed_loads[2*i-1] = DistributedLoads(assembly, i; s1=x[i], s2=x[i+1], fy = (s) -> interp_aero_force(s, x, Tp))
        distributed_loads[2*i] = DistributedLoads(assembly, i; s1=x[i], s2=x[i+1], fz = (s) -> interp_aero_force(s, x, Np))
    end

    system, converged = steady_state_analysis(assembly;
        prescribed_conditions=bcs,
        distributed_loads=distributed_loads,
        angular_velocity=[0.0, 0.0, omega], linear=true)

    return assembly, system
end

function SolverComp(; rho, E, nu, x)

    # create the points
    nnodes = length(x)
    nelems = nnodes-1
    y = zero(x)
    z = zero(x)
    points = [[x[i], y[i], z[i]] for i = 1:nnodes]

    function apply_nonlinear_forwarddiffable!(y, x)
        T = eltype(x)

        # Unpack the inputs.
        omega = x[:omega]
        Np = x[:Np]
        Tp = x[:Tp]
        chord = x[:chord]
        twist = x[:twist]
        A = x[:A]
        Iyy = x[:Iyy]
        Izz = x[:Izz]
        Iyz = x[:Iyz]

        assembly, system = create_assembly(points=self.points, omega=omega,
                                           Tp=Tp, Np=Np, chord=chord, twist=twist,
                                           A=A, Iyy=Iyy, Izz=Izz, Iyz=Iyz,
                                           rho=rho, E=E, nu=nu)

        bcs = Dict(
            1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
        )
        state = AssemblyState(system, assembly; prescribed_conditions=bcs)

        # Put the outputs in the output array.
        y[:Fx] .= [state.elements[ielem].F[1] for ielem = 1:length(assembly.elements)]
        y[:My] .= [state.elements[ielem].M[2] for ielem = 1:length(assembly.elements)]
        y[:Mz] .= [state.elements[ielem].M[3] for ielem = 1:length(assembly.elements)]

        return nothing
    end

    # Initialize the input and output vectors needed by ForwardDiff.jl. (The
    # ForwardDiff.jl inputs include phi, but that's an OpenMDAO output.)
    X = ComponentArray(
        omega=0.0, Np=zeros(Float64, nnodes), Tp=zeros(Float64, nnodes),
        chord=zeros(Float64, nelems), twist=zeros(Float64, nelems),
        A=zeros(Float64, nelems), Iyy=zeros(Float64, nelems),
        Izz=zeros(Float64, nelems), Iyz=zeros(Float64, nelems))
    Y = ComponentArray(
        Fx=zeros(Float64, nelems), My=zeros(Float64, nelems), Mz=zeros(Float64, nelems))
    J = Y.*X'

    # Get the JacobianConfig object, which we'll reuse each time when calling
    # the ForwardDiff.jacobian! function (apparently good for efficiency).
    config = ForwardDiff.JacobianConfig(apply_nonlinear_forwarddiffable!, Y, X)

    return SolverComp(rho, E, nu, points, apply_nonlinear_forwarddiffable!, X, Y, J, config)
end

function OpenMDAO.setup(self::SolverComp)

    nelems = length(self.points)-1

    # Declare the OpenMDAO inputs.
    input_data = Vector{VarData}()
    push!(input_data, VarData("omega", shape=1, units="rad/s"))
    push!(input_data, VarData("Tp", shape=nelems+1, units="N/m"))
    push!(input_data, VarData("Np", shape=nelems+1, units="N/m"))
    push!(input_data, VarData("chord", shape=nelems, units="m"))
    push!(input_data, VarData("twist", shape=nelems, units="rad"))
    push!(input_data, VarData("A", shape=nelems, units="m**2"))
    push!(input_data, VarData("Iyy", shape=nelems, units="m**4"))
    push!(input_data, VarData("Izz", shape=nelems, units="m**4"))
    push!(input_data, VarData("Iyz", shape=nelems, units="m**4"))

    # Declare the OpenMDAO outputs.
    output_data = Vector{VarData}()
    push!(output_data, VarData("Fx", shape=shape=nelems, units="N/m"))
    push!(output_data, VarData("My", shape=shape=nelems, units="N/m"))
    push!(output_data, VarData("Mz", shape=shape=nelems, units="N/m"))

    # Declare the OpenMDAO partial derivatives.
    partials_data = Vector{PartialsData}()
    push!(partials_data, PartialsData("Fx", "omega"))
    push!(partials_data, PartialsData("Fx", "A"))

    push!(partials_data, PartialsData("My", "omega"))
    push!(partials_data, PartialsData("My", "Tp"))
    push!(partials_data, PartialsData("My", "Np"))
    push!(partials_data, PartialsData("My", "A"))
    push!(partials_data, PartialsData("My", "Iyy"))
    push!(partials_data, PartialsData("My", "Izz"))
    push!(partials_data, PartialsData("My", "Iyz"))

    push!(partials_data, PartialsData("Mz", "omega"))
    push!(partials_data, PartialsData("Mz", "Tp"))
    push!(partials_data, PartialsData("Mz", "Np"))
    push!(partials_data, PartialsData("Mz", "A"))
    push!(partials_data, PartialsData("Mz", "Iyy"))
    push!(partials_data, PartialsData("Mz", "Izz"))
    push!(partials_data, PartialsData("Mz", "Iyz"))

    return input_data, output_data, partials_data
end

function OpenMDAO.compute!(self::SolverComp, inputs, outputs)

    # Unpack the inputs.
    omega = inputs["omega"][1]
    Np = inputs["Np"]
    Tp = inputs["Tp"]
    chord = inputs["chord"]
    twist = inputs["twist"]
    A = inputs["A"]
    Iyy = inputs["Iyy"]
    Izz = inputs["Izz"]
    Iyz = inputs["Iyz"]

    # Unpack the outputs.
    Fx = outputs["Fx"]
    My = outputs["My"]
    Mz = outputs["Mz"]

    E = self.E
    rho = self.rho
    nu = self.nu

    assembly, system = create_assembly(points=self.points, omega=omega,
                                       Tp=Tp, Np=Np, chord=chord, twist=twist,
                                       A=A, Iyy=Iyy, Izz=Izz, Iyz=Iyz,
                                       rho=rho, E=E, nu=nu)

    bcs = Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
    )
    state = AssemblyState(system, assembly; prescribed_conditions=bcs)

    Fx[1:end] = [state.elements[ielem].F[1] for ielem = 1:length(assembly.elements)]
    My[1:end] = [state.elements[ielem].M[2] for ielem = 1:length(assembly.elements)]
    Mz[1:end] = [state.elements[ielem].M[3] for ielem = 1:length(assembly.elements)]

end

function OpenMDAO.compute_partials!(self::SolverComp, inputs, partials)

    # Unpack the inputs.
    omega = inputs["omega"][1]
    Np = inputs["Np"]
    Tp = inputs["Tp"]
    chord = inputs["chord"]
    twist = inputs["twist"]
    A = inputs["A"]
    Iyy = inputs["Iyy"]
    Izz = inputs["Izz"]
    Iyz = inputs["Iyz"]
    nelems = length(A)
    nnodes = nelems+1

    x_ce = ComponentArray(omega=omega, Np=Np, Tp=Tp, chord=chord, twist=twist, A=A, Iyy=Iyy, Izz=Izz, Iyz=Iyz)

    # Working arrays and configuration for ForwardDiff's Jacobian routine.
    x = self.x
    y = self.y
    J = self.J
    config = self.forwarddiff_config

    x[:omega] = omega
    x[:Np] .= Np
    x[:Tp] .= Tp
    x[:chord] .= chord
    x[:twist] .= twist
    x[:A] .= A
    x[:Iyy] .= Iyy
    x[:Izz] .= Izz
    x[:Iyz] .= Iyz

    # Reshape the partials
    dFx_domega = partials["Fx", "omega"]
    dFx_dA = partials["Fx", "A"] #transpose(reshape(partials["Fx", "A"], nelems, nelems))

    dMy_domega = partials["My", "omega"]
    dMy_dTp = partials["My", "A"] #transpose(reshape(partials["My", "Tp"], nelems, nnodes))
    dMy_dNp = partials["My", "Np"] #transpose(reshape(partials["My", "Np"], nelems, nnodes))
    dMy_dA = partials["My", "A"] #transpose(reshape(partials["My", "A"], nelems, nelems))
    dMy_dIyy = partials["My", "Iyy"] #transpose(reshape(partials["My", "Iyy"], nelems, nelems))
    dMy_dIzz = partials["My", "Izz"] #transpose(reshape(partials["My", "Izz"], nelems, nelems))
    dMy_dIyz = partials["My", "Iyz"] #transpose(reshape(partials["My", "Iyz"], nelems, nelems))

    dMz_domega = partials["Mz", "omega"]
    dMz_dTp = partials["My", "Tp"] #transpose(reshape(partials["Mz", "Tp"], nelems, nnodes))
    dMz_dNp = partials["My", "Np"] #transpose(reshape(partials["Mz", "Np"], nelems, nnodes))
    dMz_dA = partials["My", "A"] #transpose(reshape(partials["Mz", "A"], nelems, nelems))
    dMz_dIyy = partials["My", "Iyy"] #transpose(reshape(partials["Mz", "Iyy"], nelems, nelems))
    dMz_dIzz = partials["My", "Izz"] #transpose(reshape(partials["Mz", "Izz"], nelems, nelems))
    dMz_dIyz = partials["My", "Iyz"] #transpose(reshape(partials["Mz", "Iyz"], nelems, nelems))

    # Get the Jacobian.
    ForwardDiff.jacobian!(J, self.apply_nonlinear_forwarddiffable!, y, x, config)

    dFx_domega .= J[:Fx, :omega]
    dFx_dA .= J[:Fx, :A]

    dMy_domega .= J[:My, :omega]
    dMy_dTp .= J[:My, :Tp]
    dMy_dNp .= J[:My, :Np]
    dMy_dA .= J[:My, :A]
    dMy_dIyy .= J[:My, :Iyy]
    dMy_dIzz .= J[:My, :Izz]
    dMy_dIyz .= J[:My, :Iyz]

    dMz_domega .= J[:Mz, :omega]
    dMz_dTp .= J[:Mz, :Tp]
    dMz_dNp .= J[:Mz, :Np]
    dMz_dA .= J[:Mz, :A]
    dMz_dIyy .= J[:Mz, :Iyy]
    dMz_dIzz .= J[:Mz, :Izz]
    dMz_dIyz .= J[:Mz, :Iyz]

end