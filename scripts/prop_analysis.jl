using DataFrames
using CSV

using GXBeam, LinearAlgebra
using CCBladeLoadingExample


function NACA0012(z)
    # For relative chord dimension z in [0, 1], return the +/- relative
    # thickness coordinates of a NACA 0012 airfoil

    y = 0.594689181*[0.298222773 *(z.^0.5) - 0.127125232 .*z - 0.357907906 .*(z.^2) + 0.291984971 .*(z.^3) - 0.105174606 .*(z.^4)]

    return y[1], -1 .*y[1]
end

function create_assembly(; points, Tp, Np, omega, chord, theta, A, Iyy, Izz, Iyz, rho, E, nu, linear=true)

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
        distributed_loads[2*i-1] = DistributedLoads(assembly, i; s1=x[i], s2=x[i+1], fy = (s) -> interp_aero_force(s, x, Np))
        distributed_loads[2*i] = DistributedLoads(assembly, i; s1=x[i], s2=x[i+1], fz = (s) -> interp_aero_force(s, x, Tp))
    end

    system, converged = steady_state_analysis(assembly;
        prescribed_conditions=bcs,
        distributed_loads=distributed_loads,
        angular_velocity=[0.0, omega, 0.0], linear=linear)

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

function run_analysis(;linear=true, chord=nothing, theta=nothing, fname="prop")

    # Compute the aerodynamic forces and set up the element discretization
    x, Np, Tp = CCBladeLoadingExample.run_ccblade()
    nelems = length(x)

    # Initialize the design variables
    span = 12.0*0.0254
    Rhub = 0.2*span
    if chord == nothing
        chord = 1.0*0.0254*ones(nelems)  # (inch)
    end
    if theta == nothing
        theta = (40.0*pi/180.0)*ones(nelems)  # (rad)
    end

    # Define the angular rotation
    rpm = 7110.0
    omega = rpm*2*pi/60

    #Tp = zero(Tp)

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
    xvals = collect(range(Rhub, span, length=nelems+1))
    for i = 1:nelems+1
        xpts[i] = xvals[i]
    end

    ypts = zero(xpts)
    zpts = zero(xpts)
    points = [[xpts[i], ypts[i], zpts[i]] for i = 1:nelems+1]

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

    #print(maximum(My))

    write_output(assembly, state; points=points, chord=chord, theta=theta, fname=fname)

    return
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

# Get saved design values to pass in to the analysis
#csv_name = "chord_theta_SNOPT_w_splines_no_ks_sf_4.0.csv"
csv_name = "coupled_chord_theta_sf_0.0.csv"
df = DataFrame(CSV.File(csv_name))
chord = df[:, :chord]
theta = df[:, :theta]

output_only(chord, theta; fname="coupled_chord_theta_sf_0")
#run_analysis(linear=true, chord=chord, theta=theta, fname="aero_opt")