using ConcreteStructs
using ComponentArrays
using ForwardDiff
using StaticArrays
using GXBeam, LinearAlgebra
using OpenMDAO: AbstractExplicitComp

@concrete struct SolverComp <: AbstractExplicitComp
    rho
    E
    nu
    Rhub
    span
    nelems
    apply_nonlinear_forwarddiffable!
    x
    y
    J
    forwarddiff_config
end

function create_assembly(; points, Tp, Np, omega, chord, theta, A, Iyy, Izz, Iyz, rho, E, nu)

    x = [points[i, 1][1] for i = 1:length(points)]
    nelems = length(points)-1
    start = 1:nelems
    stop = 2:nelems+1

    G = E/(2*(1 + nu))

    x3c = 0.17*cos.(theta).*chord
    ky = 10*(1 + nu)/(12 + 11*nu)
    kz = 10*(1 + nu)/(12 + 11*nu)

    H22 = E*Iyy
    H33 = E*Izz
    H23 = E*Iyz

    Ay = A./ky
    Az = A./kz

    T = eltype(A)
    stiffness = fill(zeros(T, (6,6)), nelems)
    mass = fill(zeros(T, (6, 6)), nelems)

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

    # Convert the stiffness matrices to compliance matrices of the same type
    assembly = Assembly(points, start, stop, stiffness=stiffness, mass=mass)

    # set prescribed conditions (fixed left endpoint)
    bcs = Dict(
        1 => PrescribedConditions(ux=zero(T), uy=zero(T), uz=zero(T), theta_x=zero(T), theta_y=zero(T), theta_z=zero(T))
    )

    function interp_aero_force(x, radii, f_aero)
        # Linearly interpolate the aerodynamic loads
        # x: coordinate where we want to evaluate the force
        # radii: x locations of nodes
        # f_aero: vector with aero forces

        if x < radii[1] || x >= radii[end]  # x not in range of radii
            q = 0
        elseif x == radii[1]
            q = f_aero[1]
        else
            k = findfirst(radii .>= x)
            j = k - 1
            q = f_aero[j]
        end

        return q
    end

    # define distributed loads
    distributed_loads = Dict()
    for i = 1:nelems
        distributed_loads[2*i-1] = DistributedLoads(assembly, i; s1=x[i], s2=x[i+1], fy = (s) -> interp_aero_force(s, x, Np))
        distributed_loads[2*i] = DistributedLoads(assembly, i; s1=x[i], s2=x[i+1], fz = (s) -> interp_aero_force(s, x, Tp))
    end

    system, converged = steady_state_analysis(assembly;
        prescribed_conditions=bcs,
        distributed_loads=distributed_loads,
        angular_velocity=[zero(T), omega, zero(T)], linear=true)

    return assembly, system
end

function SolverComp(; rho, E, nu, Rhub, span, nelems)

    function apply_nonlinear_forwarddiffable!(y, x)
        T = eltype(x)

        # Unpack the inputs.
        omega = x[:omega]
        Np = x[:Np]
        Tp = x[:Tp]
        chord = x[:chord]
        theta = x[:theta]
        A = x[:A]
        Iyy = x[:Iyy]
        Izz = x[:Izz]
        Iyz = x[:Iyz]

        xpts = zeros(T, nelems+1)
        xvals = collect(range(Rhub, span, length=nelems+1))
        for i = 1:(nelems+1)
            xpts[i] = xvals[i]
        end

        ypts = zero(xpts)
        zpts = zero(xpts)
        points = [[xpts[i], ypts[i], zpts[i]] for i = 1:(nelems+1)]

        assembly, system = create_assembly(points=points, omega=omega,
                                           Tp=Tp, Np=Np, chord=chord, theta=theta,
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

    # Initialize the input and output vectors needed by ForwardDiff.jl.
    X = ComponentArray(
        omega=0.0, Np=zeros(Float64, nelems), Tp=zeros(Float64, nelems),
        chord=zeros(Float64, nelems), theta=zeros(Float64, nelems),
        A=zeros(Float64, nelems), Iyy=zeros(Float64, nelems),
        Izz=zeros(Float64, nelems), Iyz=zeros(Float64, nelems))
    Y = ComponentArray(
        Fx=zeros(Float64, nelems), My=zeros(Float64, nelems), Mz=zeros(Float64, nelems))
    J = Y.*X'

    # Get the JacobianConfig object, which we'll reuse each time when calling
    # the ForwardDiff.jacobian! function (apparently good for efficiency).
    config = ForwardDiff.JacobianConfig(apply_nonlinear_forwarddiffable!, Y, X)

    return SolverComp(rho, E, nu, Rhub, span, nelems, apply_nonlinear_forwarddiffable!, X, Y, J, config)
end

function OpenMDAO.setup(self::SolverComp)

    nelems = self.nelems

    # Declare the OpenMDAO inputs.
    input_data = Vector{VarData}()
    push!(input_data, VarData("omega", shape=1, units="rad/s"))
    push!(input_data, VarData("Tp", shape=nelems, units="N/m"))
    push!(input_data, VarData("Np", shape=nelems, units="N/m"))
    push!(input_data, VarData("chord", shape=nelems, units="m"))
    push!(input_data, VarData("theta", shape=nelems, units="rad"))
    push!(input_data, VarData("A", shape=nelems, units="m**2"))
    push!(input_data, VarData("Iyy", shape=nelems, units="m**4"))
    push!(input_data, VarData("Izz", shape=nelems, units="m**4"))
    push!(input_data, VarData("Iyz", shape=nelems, units="m**4"))

    # Declare the OpenMDAO outputs.
    output_data = Vector{VarData}()
    push!(output_data, VarData("Fx", shape=shape=nelems, units="N"))
    push!(output_data, VarData("My", shape=shape=nelems, units="N*m"))
    push!(output_data, VarData("Mz", shape=shape=nelems, units="N*m"))

    # Declare the OpenMDAO partial derivatives.
    partials_data = Vector{PartialsData}()

    push!(partials_data, PartialsData("Fx", "omega"))
    push!(partials_data, PartialsData("Fx", "A"))

    push!(partials_data, PartialsData("My", "omega"))
    push!(partials_data, PartialsData("My", "Tp"))
    push!(partials_data, PartialsData("My", "Np"))
    push!(partials_data, PartialsData("My", "chord"))
    push!(partials_data, PartialsData("My", "theta"))
    push!(partials_data, PartialsData("My", "A"))
    push!(partials_data, PartialsData("My", "Iyy"))
    push!(partials_data, PartialsData("My", "Izz"))
    push!(partials_data, PartialsData("My", "Iyz"))

    push!(partials_data, PartialsData("Mz", "omega"))
    push!(partials_data, PartialsData("Mz", "Tp"))
    push!(partials_data, PartialsData("Mz", "Np"))
    push!(partials_data, PartialsData("Mz", "chord"))
    push!(partials_data, PartialsData("Mz", "theta"))
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
    theta = inputs["theta"]
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

    xpts = collect(range(self.Rhub, self.span, length=self.nelems+1))
    ypts = zero(xpts)
    zpts = zero(xpts)
    points = [[xpts[i], ypts[i], zpts[i]] for i = 1:(self.nelems+1)]

    assembly, system = create_assembly(points=points, omega=omega,
                                       Tp=Tp, Np=Np, chord=chord, theta=theta,
                                       A=A, Iyy=Iyy, Izz=Izz, Iyz=Iyz,
                                       rho=rho, E=E, nu=nu)

    bcs = Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
    )
    state = AssemblyState(system, assembly; prescribed_conditions=bcs)

    Fx .= [state.elements[ielem].F[1] for ielem = 1:length(assembly.elements)]
    My .= [state.elements[ielem].M[2] for ielem = 1:length(assembly.elements)]
    Mz .= [state.elements[ielem].M[3] for ielem = 1:length(assembly.elements)]

end

function OpenMDAO.compute_partials!(self::SolverComp, inputs, partials)

    nelems = self.nelems

    # Unpack the inputs.
    omega = inputs["omega"][1]
    Np = inputs["Np"]
    Tp = inputs["Tp"]
    chord = inputs["chord"]
    theta = inputs["theta"]
    A = inputs["A"]
    Iyy = inputs["Iyy"]
    Izz = inputs["Izz"]
    Iyz = inputs["Iyz"]

    # Working arrays and configuration for ForwardDiff's Jacobian routine.
    x = self.x
    y = self.y
    J = self.J
    config = self.forwarddiff_config

    x[:omega] = omega
    x[:Np] .= Np
    x[:Tp] .= Tp
    x[:chord] .= chord
    x[:theta] .= theta
    x[:A] .= A
    x[:Iyy] .= Iyy
    x[:Izz] .= Izz
    x[:Iyz] .= Iyz

    # Reshape the partials
    dFx_domega = partials["Fx", "omega"]
    dFx_dA = partials["Fx", "A"]

    dMy_domega = partials["My", "omega"]
    dMy_dTp = partials["My", "Tp"]
    dMy_dNp = partials["My", "Np"]
    dMy_dchord = partials["My", "chord"]
    dMy_dtheta = partials["My", "theta"]
    dMy_dA = partials["My", "A"]
    dMy_dIyy = partials["My", "Iyy"]
    dMy_dIzz = partials["My", "Izz"]
    dMy_dIyz = partials["My", "Iyz"]

    dMz_domega = partials["Mz", "omega"]
    dMz_dTp = partials["Mz", "Tp"]
    dMz_dNp = partials["Mz", "Np"]
    dMz_dchord = partials["Mz", "chord"]
    dMz_dtheta = partials["Mz", "theta"]
    dMz_dA = partials["Mz", "A"]
    dMz_dIyy = partials["Mz", "Iyy"]
    dMz_dIzz = partials["Mz", "Izz"]
    dMz_dIyz = partials["Mz", "Iyz"]

    # Get the Jacobian.
    ForwardDiff.jacobian!(J, self.apply_nonlinear_forwarddiffable!, y, x, config)

    dFx_domega .= J[:Fx, :omega]
    dFx_dA .= J[:Fx, :A]

    dMy_domega .= J[:My, :omega]
    dMy_dTp .= J[:My, :Tp]
    dMy_dNp .= J[:My, :Np]
    dMy_dA .= J[:My, :A]
    dMy_dchord .= J[:My, :chord]
    dMy_dtheta .= J[:My, :theta]
    dMy_dIyy .= J[:My, :Iyy]
    dMy_dIzz .= J[:My, :Izz]
    dMy_dIyz .= J[:My, :Iyz]

    dMz_domega .= J[:Mz, :omega]
    dMz_dTp .= J[:Mz, :Tp]
    dMz_dNp .= J[:Mz, :Np]
    dMz_dchord .= J[:Mz, :chord]
    dMz_dtheta .= J[:Mz, :theta]
    dMz_dA .= J[:Mz, :A]
    dMz_dIyy .= J[:Mz, :Iyy]
    dMz_dIzz .= J[:Mz, :Izz]
    dMz_dIyz .= J[:Mz, :Iyz]

end