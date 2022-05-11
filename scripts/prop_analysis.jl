using DataFrames
using CSV
using Plots
using KrylovKit
using LinearMaps

using GXBeam, LinearAlgebra
using CCBladeLoadingExample


function gram_schmidt(V, j)
    # Using the Gram-Schmidt algorithm, construct a new basis vector
    # orthogonal to each of the first j columns of V

    n = size(V, 1)

    # Make an initial guess of the vector, using a sine function for no reason
    x = collect(range(0.0, pi, length=n))
    v = sin.(x)
    u = v
    println(v)
    println(V[:, 1])
    for i in 1:j
        u -= V[:, i]*dot(V[:, i], v)/dot(V[:, i], V[:, i])
        println(u)
    end

    return u/norm(u)
end

function lanczos(A, m)
    # Lanczos algorithm to compute the m lowest eigenvalues of matrix A

    n = size(A, 1)
    T = zeros((m, m))
    V = zeros((n, m))

    # Initialize the first vector
    v = collect(range(0.0, 1.0, length=n))
    v = v/norm(v)
    V[:, 1] .= v

    # Compute the first step
    w = A*v
    a = dot(w, v)
    w = w - a*v
    T[1, 1] = a
    for j in 2:m
        # Compute b and store in in T
        b = norm(w)
        T[j, j-1] = b
        T[j-1, j] = b

        # Compute the new v
        if abs(b) > 1e-8
            v = w/b
        else
            println("using gram schmidt")
            v = gram_schmidt(V, j-1)
        end
        V[:, j] .= v

        # Compute the update to w
        w = A*v
        a = dot(w, v)
        w = w - a*v - b*V[:, j-1]

        # Add the new a to T
        T[j, j] = a
    end

    display(T)
    F = eigvals(T)

    return F
end


function arnoldi(A, m, sigma=1.0)
    # Find the m smallest eigenvalues of A using the Arnoldi algorithm

    M = inv(A + sigma*I)

    n = size(A, 1)
    b = collect(range(0.0, 1.0, length=n))
    b /= norm(b)

    eps = 1e-10
    h = zeros((m+1, m))
    Q = zeros((n, m+1))

    Q[:, 1] = b/norm(b)
    for k in 2:m+1
        v = A*Q[:, k-1]
        for j in 1:k
            h[j, k-1] = dot(Q[:, j], v)
            v = v - h[j, k-1] * Q[:, j]
        end

        h[k, k-1] = norm(v)
        if h[k, k-1] > eps
            Q[:, k] = v/h[k, k-1]
        else
            H = Q'*A*Q
            y = eigvals(H)
            return y
        end
    end

    H = Q'*A*Q
    y = eigvals(H)

    return y
end


function run_frequency_analysis(; chord=nothing, theta=nothing, omega=nothing, Np=nothing, Tp=nothing)

    # Compute the aerodynamic forces and set up the element discretization
    xe_a, Np0, Tp0 = CCBladeLoadingExample.run_ccblade()

    if (Np == nothing) & (Tp == nothing)
        nelems_aero = length(xe_a)
    elseif Np == nothing
        nelems_aero = length(Tp)
    else
        nelems_aero = length(Np)
    end

    if Np == nothing
        Np = Np0
    end
    if Tp == nothing
        Tp = Tp0
    end

    if (chord == nothing) & (theta == nothing)
        nelems = length(xe_a)
    elseif chord == nothing
        nelems = length(theta)
    else
        nelems = length(chord)
    end

    span = 12.0*0.0254
    Rhub = 0.2*span
    xvals = collect(range(Rhub, span, length=nelems+1))
    xvals_aero = collect(range(Rhub, span, length=nelems_aero+1))

    # Initialize the design variables
    if chord == nothing
        chord = 1.0*0.0254*ones(nelems)  # (inch)
    end
    if theta == nothing
        theta = (40.0*pi/180.0)*ones(nelems)  # (rad)
    end

    # Define the angular rotation
    if omega == nothing
        rpm = 7200.0
        omega = rpm*2*pi/60
    end

    # Define the reference area properties
    A_ref = 821.8
    Iyy_ref = 23543.4
    Izz_ref = 5100.8

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

    xpts = zeros(nelems+1)
    for i = 1:nelems+1
        xpts[i] = xvals[i]
    end

    ypts = zero(xpts)
    zpts = zero(xpts)
    points = [[xpts[i], ypts[i], zpts[i]] for i = 1:nelems+1]

    assembly, system = create_assembly(points=points, x_aero=xvals_aero, omega=omega,
                                       Tp=Tp, Np=Np, chord=chord, theta=theta,
                                       A=A, Iyy=Iyy, Izz=Izz, Iyz=Iyz,
                                       rho=rho, E=E, nu=nu)

    # set prescribed conditions (fixed left endpoint)
    bcs = Dict(
        1 => PrescribedConditions(ux=0.0, uy=0.0, uz=0.0, theta_x=0.0, theta_y=0.0, theta_z=0.0)
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
        distributed_loads[2*i-1] = DistributedLoads(assembly, i; s1=xvals_aero[i], s2=xvals_aero[i+1], fy = (s) -> interp_aero_force(s, xpts, Np))
        distributed_loads[2*i] = DistributedLoads(assembly, i; s1=xvals_aero[i], s2=xvals_aero[i+1], fz = (s) -> interp_aero_force(s, xpts, Tp))
    end

    # system, λ, V, converged = eigenvalue_analysis!(system, assembly;
    #         prescribed_conditions=bcs,
    #         distributed_loads=distributed_loads,
    #         angular_velocity=[0.0, omega, 0.0],
    #         nev=6)
    # println(imag(λ))
    # println("eigenvalue analysis converged: ", converged)

    K = system.K
    M = system.M
    x = system.x

    # unpack scaling parameter
    force_scaling = system.force_scaling

    # also unpack system indices
    irow_point = system.irow_point
    irow_elem = system.irow_elem
    irow_elem1 = system.irow_elem1
    irow_elem2 = system.irow_elem2
    icol_point = system.icol_point
    icol_elem = system.icol_elem

    # current time
    t = system.t

    # current parameters
    pcond = bcs
    dload = distributed_loads
    pmass = Dict{Int,PointMass{Float64}}()  #typeof(point_masses) <: AbstractDict ? point_masses : point_masses(t)
    gvec = [0.0, 0.0, 0.0]  #typeof(gravity) <: AbstractVector ? SVector{3}(gravity) : SVector{3}(gravity(tvec[it]))
    x0 = [0.0, 0.0, 0.0]  #typeof(origin) <: AbstractVector ? SVector{3}(origin) : SVector{3}(origin(t))
    v0 = [0.0, 0.0, 0.0]  #typeof(linear_velocity) <: AbstractVector ? SVector{3}(linear_velocity) : SVector{3}(linear_velocity(t))
    ω0 = [0.0, omega, 0.0]  #typeof(angular_velocity) <: AbstractVector ? SVector{3}(angular_velocity) : SVector{3}(angular_velocity(t))
    a0 = [0.0, 0.0, 0.0]  #typeof(linear_acceleration) <: AbstractVector ? SVector{3}(linear_acceleration) : SVector{3}(linear_acceleration(t))
    α0 = [0.0, 0.0, 0.0]  #typeof(angular_acceleration) <: AbstractVector ? SVector{3}(angular_acceleration) : SVector{3}(angular_acceleration(t))

    # solve for the system stiffness matrix
    K = GXBeam.steady_state_system_jacobian!(K, x, assembly, pcond, dload, pmass, gvec, force_scaling,
        irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0,
        ω0, a0, α0)

    # solve for the system mass matrix
    M = GXBeam.system_mass_matrix!(M, x, assembly, pmass, force_scaling, irow_point,
        irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem)

    A = K*inv(M)

#     # construct linear map
#     T = eltype(system)
#     nx = length(x)
#     Kfact = GXBeam.safe_lu(K)
#     f! = (b, x) -> ldiv!(b, Kfact, M * x)
#     fc! = (b, x) -> mul!(b, M', Kfact' \ x)
#     A = LinearMap{T}(f!, fc!, nx, nx; ismutating=true)
#
#     # # Solve the generalized eigenvalue problem
#     # vals, vecs, info = KrylovKit.geneigsolve((K, M), howmany=6)
#     vals, vecs, info = KrylovKit.eigsolve(A, x, howmany=6, which=:LM)
#     # println(vals)
#     # println(info)

end


function compare_linear_nonlinear(xe; chord=nothing, theta=nothing, Np=nothing, Tp=nothing)

    state, assembly = run_analysis(linear=true, chord=chord, theta=theta, Np=Np, Tp=Tp)
    Fx = [state.elements[ielem].F[1] for ielem = 1:length(assembly.elements)]
    My = [state.elements[ielem].M[2] for ielem = 1:length(assembly.elements)]
    Mz = [state.elements[ielem].M[3] for ielem = 1:length(assembly.elements)]
    ux = [state.elements[ielem].u[1] for ielem = 1:length(assembly.elements)]./0.0254
    uy = [state.elements[ielem].u[2] for ielem = 1:length(assembly.elements)]./0.0254
    uz = [state.elements[ielem].u[3] for ielem = 1:length(assembly.elements)]./0.0254

    state_nl, assembly = run_analysis(linear=false, chord=chord, theta=theta, Np=Np, Tp=Tp)
    Fx_nl = [state_nl.elements[ielem].F[1] for ielem = 1:length(assembly.elements)]
    My_nl = [state_nl.elements[ielem].M[2] for ielem = 1:length(assembly.elements)]
    Mz_nl = [state_nl.elements[ielem].M[3] for ielem = 1:length(assembly.elements)]
    ux_nl = [state_nl.elements[ielem].u[1] for ielem = 1:length(assembly.elements)]./0.0254
    uy_nl = [state_nl.elements[ielem].u[2] for ielem = 1:length(assembly.elements)]./0.0254
    uz_nl = [state_nl.elements[ielem].u[3] for ielem = 1:length(assembly.elements)]./0.0254

    l = @layout [a ; b; c]
    p1 = plot(xe, Fx,
              xlim = (xe[1], xe[length(xe)]),
              xticks = 0.0:1.0:xe[length(xe)],
              xlabel = "x (in)",
              ylabel = "Fx (N)",
              grid = false,
              label="Linear")
    plot!(xe, Fx_nl, label="Nonlinear")
    p2 = plot(xe, My,
              xlim = (xe[1], xe[length(xe)]),
              xticks = 0.0:1.0:xe[length(xe)],
              xlabel = "x (in)",
              ylabel = "My (N-m)",
              grid = false,
              label="Linear")
    plot!(xe, My_nl, label="Nonlinear")
    p3 = plot(xe, Mz,
              xlim = (xe[1], xe[length(xe)]),
              xticks = 0.0:1.0:xe[length(xe)],
              xlabel = "x (in)",
              ylabel = "Mz (N-m)",
              grid = false,
              label="Linear")
    plot!(xe, Mz_nl, label="Nonlinear")
    plot(p1, p2, p3, layout=l, legend=true)

    savefig("nonlinear_force_compare.pdf")

    l = @layout [a ; b; c]
    p1 = plot(xe, ux,
              xlim = (xe[1], xe[length(xe)]),
              xticks = 0.0:1.0:xe[length(xe)],
              xlabel = "x (in)",
              ylabel = "ux (in)",
              grid = false,
              label="Linear")
    plot!(xe, ux_nl, label="Nonlinear")
    p2 = plot(xe, uy,
              xlim = (xe[1], xe[length(xe)]),
              xticks = 0.0:1.0:xe[length(xe)],
              xlabel = "x (in)",
              ylabel = "uy (in)",
              grid = false,
              label="Linear")
    plot!(xe, uy_nl, label="Nonlinear")
    p3 = plot(xe, uz,
              xlim = (xe[1], xe[length(xe)]),
              xticks = 0.0:1.0:xe[length(xe)],
              xlabel = "x (in)",
              ylabel = "uz (in)",
              grid = false,
              label="Linear")
    plot!(xe, uz_nl, label="Nonlinear")
    plot(p1, p2, p3, layout=l, legend=true)

    savefig("nonlinear_disp_compare.pdf")

    return
end

function check_nonlinear(x, ux, uy, uz, chord, theta, Np, Tp; fname="nonlinear_disp_compare.pdf")

    state_nl, assembly = run_analysis(linear=false, chord=chord, theta=theta, Np=Np, Tp=Tp)
    ux_nl = [state_nl.points[inode].u[1] for inode = 1:length(assembly.points)]./0.0254
    uy_nl = [state_nl.points[inode].u[2] for inode = 1:length(assembly.points)]./0.0254
    uz_nl = [state_nl.points[inode].u[3] for inode = 1:length(assembly.points)]./0.0254

    l = @layout [a ; b; c]
    p1 = plot(x, ux,
              xlim = (x[1], x[length(x)]),
              xticks = x[1]:1.0:x[length(x)],
              xlabel = "x (in)",
              ylabel = "ux (in)",
              grid = false,
              label="Optimization")
    plot!(x, ux_nl, label="Nonlinear analysis")
    p2 = plot(x, uy,
              xlim = (x[1], x[length(x)]),
              xticks = x[1]:1.0:x[length(x)],
              xlabel = "x (in)",
              ylabel = "uy (in)",
              grid = false,
              label="Optimization")
    plot!(x, uy_nl, label="Nonlinear analysis")
    p3 = plot(x, uz,
              xlim = (x[1], x[length(x)]),
              xticks = x[1]:1.0:x[length(x)],
              xlabel = "x (in)",
              ylabel = "uz (in)",
              grid = false,
              label="Optimization")
    plot!(x, uz_nl, label="Nonlinear analysis")
    plot(p1, p2, p3, layout=l, legend=true)

    savefig(fname)

    return
end

function NACA0012(z)
    # For relative chord dimension z in [0, 1], return the +/- relative
    # thickness coordinates of a NACA 0012 airfoil

    y = 0.594689181*[0.298222773 *(z.^0.5) - 0.127125232 .*z - 0.357907906 .*(z.^2) + 0.291984971 .*(z.^3) - 0.105174606 .*(z.^4)]

    return y[1], -1 .*y[1]
end

function create_assembly(; points, x_aero, Tp, Np, omega, chord, theta, A, Iyy, Izz, Iyz, rho, E, nu, linear=true)

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
        1 => PrescribedConditions(ux=0.0, uy=0.0, uz=0.0, theta_x=0.0, theta_y=0.0, theta_z=0.0)
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
        distributed_loads[2*i-1] = DistributedLoads(assembly, i; s1=x_aero[i], s2=x_aero[i+1], fy = (s) -> interp_aero_force(s, x, Np))
        distributed_loads[2*i] = DistributedLoads(assembly, i; s1=x_aero[i], s2=x_aero[i+1], fz = (s) -> interp_aero_force(s, x, Tp))
    end

    system, converged = steady_state_analysis(assembly;
        prescribed_conditions=bcs,
        distributed_loads=distributed_loads,
        angular_velocity=[0.0, omega, 0.0], linear=linear)
    println("converged: ", converged)

    return assembly, system
end

function output_only(chord, theta; fname="prop")

    nelems = length(chord)
    theta = theta*(pi/180.0)

    # Define the element spacing
    span = 12.0
    Rhub = 0.2*span
    xpts = collect(range(Rhub, span, length=nelems+1))
    ypts = zero(xpts)
    zpts = zero(xpts)
    points = [[xpts[i], ypts[i], zpts[i]] for i = 1:nelems+1]

    start = 1:nelems
    stop = 2:nelems+1

    compliance = fill(Diagonal([0, 0, 0, 1.0, 1.0, 0]), nelems)

    # create assembly
    assembly = Assembly(points, start, stop, compliance=compliance)

    bcs = Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
    )

    system, converged = static_analysis(assembly;
        prescribed_conditions=bcs,
        linear=true)
    state = AssemblyState(system, assembly; prescribed_conditions=bcs)

    write_output(assembly, state; points=points, chord=chord, theta=theta, fname=fname, scaling=0.0)

    return
end

function run_analysis(;linear=true, chord=nothing, theta=nothing, omega=nothing, Np=nothing, Tp=nothing, fname="prop")

    # Compute the aerodynamic forces and set up the element discretization
    xe_a, Np0, Tp0 = CCBladeLoadingExample.run_ccblade()

    if (Np == nothing) & (Tp == nothing)
        nelems_aero = length(xe_a)
    elseif Np == nothing
        nelems_aero = length(Tp)
    else
        nelems_aero = length(Np)
    end

    if Np == nothing
        Np = Np0
    end
    if Tp == nothing
        Tp = Tp0
    end

    if (chord == nothing) & (theta == nothing)
        nelems = length(xe_a)
    elseif chord == nothing
        nelems = length(theta)
    else
        nelems = length(chord)
    end

    span = 12.0*0.0254
    Rhub = 0.2*span
    xvals = collect(range(Rhub, span, length=nelems+1))
    xvals_aero = collect(range(Rhub, span, length=nelems_aero+1))

    # Initialize the design variables
    if chord == nothing
        chord = 1.0*0.0254*ones(nelems)  # (inch)
    end
    if theta == nothing
        theta = (40.0*pi/180.0)*ones(nelems)  # (rad)
    end

    # Define the angular rotation
    if omega == nothing
        rpm = 7200.0
        omega = rpm*2*pi/60
    end

    # Define the reference area properties
    A_ref = 821.8
    Iyy_ref = 23543.4
    Izz_ref = 5100.8

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

    xpts = zeros(nelems+1)
    for i = 1:nelems+1
        xpts[i] = xvals[i]
    end

    ypts = zero(xpts)
    zpts = zero(xpts)
    points = [[xpts[i], ypts[i], zpts[i]] for i = 1:nelems+1]

    assembly, system = create_assembly(points=points, x_aero=xvals_aero, omega=omega,
                                       Tp=Tp, Np=Np, chord=chord, theta=theta,
                                       A=A, Iyy=Iyy, Izz=Izz, Iyz=Iyz,
                                       rho=rho, E=E, nu=nu, linear=linear)

    bcs = Dict(
        1 => PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)
    )
    state = AssemblyState(system, assembly; prescribed_conditions=bcs)

    Fx = [state.elements[ielem].F[1] for ielem = 1:length(assembly.elements)]
    My = [state.elements[ielem].M[2] for ielem = 1:length(assembly.elements)]
    Mz = [state.elements[ielem].M[3] for ielem = 1:length(assembly.elements)]

    #print(maximum(My))

    write_output(assembly, state; points=points, chord=chord, theta=theta, fname=fname)

    return state, assembly
end

function write_output(assembly, state; points, chord, theta, fname="prop", scaling=5.0)
    # Create a vtk file output from the GXBeam analysis

    num_airfoil_eval_pts = 30

    function generate_airfoil(chord, theta, num_airfoil_eval_pts)
        # Function to generate a vector defining a NACA 0012 airfoil
        # section with scalar chord and theta values

        airfoil = zeros(2*num_airfoil_eval_pts, 2)
        s = sin(theta)
        c = cos(theta)
        z_rel = collect(range(0.0, 1.0, length=num_airfoil_eval_pts)) .- 0.25
        y1_rel, y2_rel = NACA0012(z_rel .+ 0.25)
        for i = 1:num_airfoil_eval_pts
            z1 = chord*( z_rel[i]*c + y1_rel[i]*s)
            y1 = chord*(-z_rel[i]*s + y1_rel[i]*c)
            z2 = chord*( z_rel[i]*c + y2_rel[i]*s)
            y2 = chord*(-z_rel[i]*s + y2_rel[i]*c)
            airfoil[2*i - 1, 1] = y1
            airfoil[2*i - 1, 2] = z1
            airfoil[2*i, 1]     = y2
            airfoil[2*i, 2]     = z2
        end

        return airfoil
    end

    # Create the sections
    sections = zeros(3, 2*num_airfoil_eval_pts, length(points))
    for ip = 1:length(points)
        if ip == 1
            chord_i = chord[1]
            theta_i = theta[1]
        elseif ip == length(points)
            chord_i = chord[end]
            theta_i = theta[end]
        else
            chord_i = (chord[ip-1] + chord[ip])/2.0
            theta_i = (theta[ip-1] + theta[ip])/2.0
        end

    airfoil = generate_airfoil(chord_i, theta_i, num_airfoil_eval_pts)
    sections[1, :, ip] .= 0
    sections[2, :, ip] .= airfoil[:,1]
    sections[3, :, ip] .= airfoil[:,2]
    end

    # Write out the vtk file
    write_vtk(fname, assembly, state; sections=sections, scaling=scaling)
    #write_vtk(fname, assembly, state; scaling=1.0)

    return nothing
end

# # Test Gram-Schmidt algorithm
# V = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
# println(V)
# v = gram_schmidt(V, 3)
# println(v)
# println(dot(v, V[:, 1]))
# println(dot(v, V[:, 2]))
# println(dot(v, V[:, 3]))

# # Test eigenvalue solver
# A = [5.0 1.0 0.0 0.0 0.0;
#      0.0 3.0 1.0 0.0 -1.0;
#      1.0 0.0 2.0 1.0 0.0;
#      1.0 1.0 0.0 4.0 0.0;
#      0.0 0.0 0.0 0.0 1.0]
# F = eigvals(A)
# println(F)
# y = arnoldi(A, 2)
# println(y)

# # Get saved design values to pass in to the analysis
# csv_name = "nl_dvs_using_theta_sf_4.0.csv"
# df = DataFrame(CSV.File(csv_name))
# chord = 0.0254*df[:, :chord]
# theta = df[:, :theta]
#
# # Get saved aero forces to pass in to the analysis
# csv_name = "nl_aero_forces_using_theta_sf_4.0.csv"
# df = DataFrame(CSV.File(csv_name))
# Np = df[:, :Np]
# Tp = df[:, :Tp]
#
# # Get saved displacements to check
# csv_name = "nl_disp_using_theta_sf_4.0.csv"
# df = DataFrame(CSV.File(csv_name))
# x = df[:, :x0]
# u1 = df[:, :u1]
# u2 = df[:, :u2]
# u3 = df[:, :u3]
#
# #output_only(chord, theta; fname="coupled_chord_theta_sf_0")
# #run_analysis(linear=true)
# check_nonlinear(x, u1, u2, u3, chord, theta, Np, Tp)
#
# # Define the element spacing
# nelems = length(chord)
# span = 12.0
# Rhub = 0.2*span
# dx = (span-Rhub)/nelems
# xe = collect(range(Rhub+0.5*dx, span-0.5*dx, length=nelems))
#
# #compare_linear_nonlinear(xe; chord=chord, theta=theta, Np=Np, Tp=Tp)
# #run_frequency_analysis(chord=chord, theta=theta, Np=Np, Tp=Tp)
#
