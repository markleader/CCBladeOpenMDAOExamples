using ConcreteStructs
using ComponentArrays
using ForwardDiff
using StaticArrays
using GXBeam, LinearAlgebra
using CCBladeLoadingExample
using BenchmarkTools


function NACA0012(z)
    # For relative chord dimension z in [0, 1], return the +/- relative
    # thickness coordinates of a NACA 0012 airfoil

    y = 0.594689181*[0.298222773 *(z.^0.5) - 0.127125232 .*z - 0.357907906 .*(z.^2) + 0.291984971 .*(z.^3) - 0.105174606 .*(z.^4)]

    return y[1], -1 .*y[1]
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
        distributed_loads[2*i-1] = DistributedLoads(assembly, i; s1=x[i], s2=x[i+1], fy = (s) -> interp_aero_force(s, x, Np))
        distributed_loads[2*i] = DistributedLoads(assembly, i; s1=x[i], s2=x[i+1], fz = (s) -> interp_aero_force(s, x, Tp))
    end

    system, converged = steady_state_analysis(assembly;
        prescribed_conditions=bcs,
        distributed_loads=distributed_loads,
        angular_velocity=[zero(T), omega, zero(T)], linear=true)

    return assembly, system
end

function apply_nonlinear_forwarddiffable!(y, x)

    T = eltype(x)

    # Define the reference area properties
    A_ref = 821.8
    Iyy_ref = 23543.4
    Izz_ref = 5100.8

    # Unpack the inputs.
    omega = x[:omega]
    Np = x[:Np]
    Tp = x[:Tp]
    chord = x[:chord]
    theta = x[:theta]

    span = 12.0*0.0254
    nelems = length(chord)
    nnodes = nelems+1
    num_stress_eval_points = 20

    # Compute the area properties
    s = sin.(theta)
    c = cos.(theta)
    s2 = sin.(theta).^2
    c2 = cos.(theta).^2

    k = chord/100  # scale factor
    A = (k.^2)*A_ref
    Iyy0 = (k.^4)*Iyy_ref
    Izz0 = (k.^4)*Izz_ref
    Iyy = @. c2 *Iyy0 + s2 *Izz0
    Izz = @. s2 *Iyy0 + c2 *Izz0
    Iyz = @. (Izz0 - Iyy0)*s*c

    # Solve the finite element problem
    rho = 2780.0
    E = 72.4e9
    nu = 0.33
    ys = 345e6

    xpts = zeros(T, nnodes)
    xvals = collect(range(0.06, span, length=nnodes))
    for i = 1:nnodes
        xpts[i] = xvals[i]
    end

    ypts = zero(xpts)
    zpts = zero(xpts)
    points = [[xpts[i], ypts[i], zpts[i]] for i = 1:nnodes]

    assembly, system = create_assembly(points=points, omega=omega,
                                       Tp=Tp, Np=Np, chord=chord, theta=theta,
                                       A=A, Iyy=Iyy, Izz=Izz, Iyz=Iyz,
                                       rho=rho, E=E, nu=nu)

    bcs = Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
    )
    state = AssemblyState(system, assembly; prescribed_conditions=bcs)

    Fx = [state.elements[ielem].F[1] for ielem = 1:length(assembly.elements)]
    My = [state.elements[ielem].M[2] for ielem = 1:length(assembly.elements)]
    Mz = [state.elements[ielem].M[3] for ielem = 1:length(assembly.elements)]

    # Compute the stress in the beam
    sigma1 = zeros(T, 2*nelems*num_stress_eval_points)
    s = sin.(theta)
    c = cos.(theta)
    dI = @. Iyy*Izz - Iyz*Iyz
    z_rel = collect(range(0.0, 1.0, length=num_stress_eval_points)) .- 0.421
    y1_rel, y2_rel = NACA0012(z_rel .+ 0.421)

    for i = 1:length(chord)
        for j = 1:num_stress_eval_points
            z1 = chord[i]*(z_rel[j]*c[i] + y1_rel[j]*s[i])
            y1 = chord[i]*(-z_rel[j]*s[i] + y1_rel[j]*c[i])
            z2 = chord[i]*(z_rel[j]*c[i] + y2_rel[j]*s[i])
            y2 = chord[i]*(-z_rel[j]*s[i] + y2_rel[j]*c[i])

            idx = 2*(i-1)*num_stress_eval_points + 2*j - 1
            sigma1[idx]   = Fx[i]/A[i] - My[i]*(y1*Iyz[i] - z1*Izz[i])/dI[i] - Mz[i]*(y1*Iyy[i] - z1*Iyz[i])/dI[i]
            sigma1[idx+1] = Fx[i]/A[i] - My[i]*(y2*Iyz[i] - z2*Izz[i])/dI[i] - Mz[i]*(y2*Iyy[i] - z2*Iyz[i])/dI[i]
        end
    end
    sigma1 ./= ys

    return nothing
end


function check_totals()

    # Compute the aerodynamic forces and set up the element discretization
    x, Np, Tp = CCBladeLoadingExample.run_ccblade()
    nnodes = length(x)
    nelems = nnodes-1
    num_stress_eval_points = 20

    # Initialize the design variables
    span = 12.0*0.0254
    chord = 1.0*0.0254*ones(nelems)  # (inch)
    theta = (45.0*pi/180.0)*ones(nelems)  # (rad)

    # Define the angular rotation
    rpm = 7110.0
    omega = rpm*2*pi/60

    # Initialize the input and output vectors needed by ForwardDiff.jl.
    X = ComponentArray(
        omega=0.0, Np=zeros(Float64, nnodes), Tp=zeros(Float64, nnodes),
        chord=zeros(Float64, nelems), theta=zeros(Float64, nelems))
    Y = ComponentArray(
        sigma1=zeros(Float64, 2*nelems*num_stress_eval_points))
    J = Y.*X'

    # Get the JacobianConfig object, which we'll reuse each time when calling
    # the ForwardDiff.jacobian! function (apparently good for efficiency).
    config = ForwardDiff.JacobianConfig(apply_nonlinear_forwarddiffable!, Y, X)

    # Get the Jacobian.
    ForwardDiff.jacobian!(J, apply_nonlinear_forwarddiffable!, Y, X, config)

    return nothing
end

b = @benchmarkable check_totals() seconds=28800 samples=1
run(b)