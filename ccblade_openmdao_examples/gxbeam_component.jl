using ConcreteStructs
using GXBeam, LinearAlgebra
using OpenMDAO: AbstractExplicitComp
using Plots

@concrete struct SolverComp <: AbstractExplicitComp
    rho
    E
    nu
    points
end

function create_assembly(; points, Tp, Np, omega, chord, twist, A, Iyy, Izz, Iyz, rho, E, nu)

    x = [points[i, 1][1] for i = 1:length(points)]
    nelem = length(points)-1
    start = 1:nelem
    stop = 2:nelem+1

    G = E/(2*(1 + nu))

    x3c = 0.0 *cos.(twist).*chord  # 0.17*cos.(twist).*chord
    ky = 10*(1 + nu)/(12 + 11*nu)
    kz = 10*(1 + nu)/(12 + 11*nu)

    H22 = E*Iyy
    H33 = E*Izz
    H23 = E*Iyz

    Ay = A./ky
    Az = A./kz

    stiffness = fill(zeros(6,6), nelem)
    mass = fill(zeros(6, 6), nelem)

    for i in 1:nelem
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
    for i = 1:nelem
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
    y = zero(x)
    z = zero(x)
    points = [[x[i], y[i], z[i]] for i = 1:length(x)]

    return SolverComp(rho, E, nu, points)
end

function OpenMDAO.setup(self::SolverComp)

    nelem = length(self.points)-1

    # Declare the OpenMDAO inputs.
    input_data = Vector{VarData}()
    push!(input_data, VarData("omega", shape=1, units="rad/s"))
    push!(input_data, VarData("Tp", shape=nelem+1, units="N/m"))
    push!(input_data, VarData("Np", shape=nelem+1, units="N/m"))
    push!(input_data, VarData("chord", shape=nelem, units="m"))
    push!(input_data, VarData("twist", shape=nelem, units="rad"))
    push!(input_data, VarData("A", shape=nelem, units="m**2"))
    push!(input_data, VarData("Iyy", shape=nelem, units="m**4"))
    push!(input_data, VarData("Izz", shape=nelem, units="m**4"))
    push!(input_data, VarData("Iyz", shape=nelem, units="m**4"))

    # Declare the OpenMDAO outputs.
    output_data = Vector{VarData}()
    push!(output_data, VarData("Fx", shape=shape=nelem, units="N/m"))
    push!(output_data, VarData("My", shape=shape=nelem, units="N/m"))
    push!(output_data, VarData("Mz", shape=shape=nelem, units="N/m"))

    # Declare the OpenMDAO partial derivatives.
    partials_data = Vector{PartialsData}()
    push!(partials_data, PartialsData("Fx", "omega"))
    push!(partials_data, PartialsData("Fx", "A"))
    # push!(partials_data, PartialsData("Fx", "Iyy"))
    # push!(partials_data, PartialsData("Fx", "Izz"))
    # push!(partials_data, PartialsData("Fx", "Iyz"))
    # push!(partials_data, PartialsData("My", "A"))
    # push!(partials_data, PartialsData("My", "Iyy"))
    # push!(partials_data, PartialsData("My", "Izz"))
    # push!(partials_data, PartialsData("My", "Iyz"))
    # push!(partials_data, PartialsData("My", "omega"))
    push!(partials_data, PartialsData("My", "Tp"))
    push!(partials_data, PartialsData("My", "Np"))
    push!(partials_data, PartialsData("Mz", "Tp"))
    push!(partials_data, PartialsData("Mz", "Np"))
    # push!(partials_data, PartialsData("Mz", "A"))
    # push!(partials_data, PartialsData("Mz", "Iyy"))
    # push!(partials_data, PartialsData("Mz", "Izz"))
    # push!(partials_data, PartialsData("Mz", "Iyz"))
    # push!(partials_data, PartialsData("Mz", "omega"))

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

    dFx_domega = partials["Fx", "omega"]
    dFx_dA = partials["Fx", "A"]
    dMy_dTp = partials["My", "Tp"]
    dMy_dNp = partials["My", "Np"]
    dMz_dTp = partials["Mz", "Tp"]
    dMz_dNp = partials["Mz", "Np"]
    # dMy_dc = partials["My", "chord"]
    # dMz_dc = partials["Mz", "chord"]

    nelem = length(self.points)-1

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
    Fx = [state.elements[ielem].F[1] for ielem = 1:length(assembly.elements)]
    My = [state.elements[ielem].M[2] for ielem = 1:length(assembly.elements)]
    Mz = [state.elements[ielem].M[3] for ielem = 1:length(assembly.elements)]

    x = [self.points[i, 1][1] for i = 1:length(self.points)]
    nelem = length(self.points)-1
    start = 1:nelem
    stop = 2:nelem+1
    dl = x[stop] - x[start]
    xe = [assembly.elements[ielem].x[1] for ielem = 1:length(assembly.elements)]

    # Check that I am computing the axial forces correctly
    Fx_check = 0.5 .* (dl ./1.0) .* (omega^2 * rho .*A .* xe)

    for i = 1:length(Fx_check)
        for j = i+1:length(Fx_check)
            Fx_check[i] += 2.0 * Fx_check[j]
        end
    end
    println(abs.(Fx_check - Fx)./Fx)

    p = plot(xlim = (x[1], x[length(x)]),
        xticks = 0.0:0.05:x[length(x)],
        xlabel = "x (m)",
        ylabel = "Fx (N)",
        overwrite=true,
        grid=false,
        overwrite_figure=false,
        show=true)
    plot!(xe, Fx, label="Fx")
    plot!(xe, Fx_check, label="check")
    savefig("check_forces.pdf")

    # Compute the gradients of the axial force
    dFx_domega[1:end] = 2.0*Fx/omega

    k = rho*omega^2
    for i = 1:nelem
        for j = i:nelem
            if i == j
                dFx_dA[j, i] = 0.5*dl[i]*k*xe[i]
            else
                dFx_dA[j, i] = dl[i]*k*xe[i]
            end
        end
    end

    x = system.x
    F = system.r
    J = system.K
#     force_scaling = system.force_scaling
#     mass_scaling = system.mass_scaling
    irow_pt = system.irow_pt
    irow_beam = system.irow_beam
    irow_beam1 = system.irow_beam1
    irow_beam2 = system.irow_beam2
    icol_pt = system.icol_pt
    icol_beam = system.icol_beam
#
#     unscaled_Fx = x[icol_beam.+6]
#     force_scaling = (Fx./unscaled_Fx)[1]
#     mass_scaling = 1.0

#
#     j! = (J, x) -> system_jacobian!(J, x, assembly, prescribed_conditions,
#             distributed_loads, force_scaling, mass_scaling, irow_pt, irow_beam, irow_beam1,
#             irow_beam2, icol_pt, icol_beam, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, omega[1]])
#     j!(J, x)
#
    for i in 1:nelem
        for j in 1:nelem
            dMy_dTp[i, j] = J[irow_beam[i]+1, icol_beam[j]+10]
            dMy_dNp[i, j] = J[irow_beam[i]+2, icol_beam[j]+10]
            dMz_dTp[i, j] = J[irow_beam[i]+1, icol_beam[j]+11]
            dMz_dNp[i, j] = J[irow_beam[i]+2, icol_beam[j]+11]
        end
    end
    println(dMy_dTp)
    println(dMy_dNp)
    println(dMz_dTp)
    println(dMz_dNp)

end