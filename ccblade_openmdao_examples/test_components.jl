using OpenMDAO: om, make_component

include("area_component.jl")

prob = om.Problem()

nelems = 5
chord = collect(range(1.5*0.0254, 0.5*0.0254, length=nelems))
theta = 1.57079*ones(nelems) #collect(range(1.0, 0.0, length=nelems))

A_ref = 821.8
Iyy_ref = 23543.4
Izz_ref = 5100.8

ivc = om.IndepVarComp()
ivc.add_output("chord", chord, units="m")
ivc.add_output("theta", theta, units="rad")
prob.model.add_subsystem("ivc", ivc, promotes=["*"])

area_comp = make_component(AreaComp(nelems=nelems, A_ref=A_ref, Iyy_ref=Iyy_ref, Izz_ref=Izz_ref))  # Need to convert Julia obj to Python obj
prob.model.add_subsystem("area_comp", area_comp, promotes=["*"])

# prob.setup()
# prob.run_model()
# prob.check_partials(compact_print=true)

include("stress_component.jl")

prob = om.Problem()

nelems = 5
Nx = collect(range(1000.0, 0.0, length=nelems))
My = collect(range(100.0, 0.0, length=nelems))
Mz = collect(range(500.0, 0.0, length=nelems))
chord = collect(range(1.5*0.0254, 0.5*0.0254, length=nelems))
theta = collect(range(1.0, 0.0, length=nelems))
A = collect(range(0.01, 0.005, length=nelems))
Iyy = collect(range(0.002, 0.001, length=nelems))
Izz = collect(range(0.002, 0.001, length=nelems))
Iyz = collect(range(0.001, 0.0, length=nelems))

ivc = om.IndepVarComp()
ivc.add_output("Nx", Nx, units="N")
ivc.add_output("My", My, units="N*m")
ivc.add_output("Mz", Mz, units="N*m")
ivc.add_output("chord", chord, units="m")
ivc.add_output("theta", theta, units="rad")
ivc.add_output("A", A, units="m**2")
ivc.add_output("Iyy", Iyy, units="m**4")
ivc.add_output("Izz", Izz, units="m**4")
ivc.add_output("Iyz", Iyz, units="m**4")
prob.model.add_subsystem("ivc", ivc, promotes=["*"])

stress_comp = make_component(StressComp(nelems=nelems, num_stress_eval_points=50, ys=345e6))
prob.model.add_subsystem("stress_comp", stress_comp, promotes=["*"])

# prob.setup(force_alloc_complex=true)
# prob.run_model()
# prob.check_partials(compact_print=true, method="cs")

include("gxbeam_component.jl")

prob = om.Problem()

span = 0.3
Rhub = 0.2*span
nelems = 10
omega = 7110.0*2*pi/60
xe = collect(range(Rhub, span, length=nelems))
Tp = collect(range(10.0, 50.0, length=nelems))
Np = collect(range(10.0, 200.0, length=nelems))
chord = collect(range(1.5*0.0254, 0.5*0.0254, length=nelems))
theta = collect(range(1.0, 0.0, length=nelems))

s = sin.(theta)
c = cos.(theta)
s2 = s.^2
c2 = c.^2

k = chord/100  # scale factor
A = (k.^2)*821.8
Iyy0 = (k.^4)*23543.4
Izz0 = (k.^4)*5100.8
Iyy = @. c2 *Iyy0 + s2 *Izz0
Izz = @. s2 *Iyy0 + c2 *Izz0
Iyz = @. (Iyy0 - Izz0)*s*c

ivc = om.IndepVarComp()
ivc.add_output("omega", omega, units="rad/s")
ivc.add_output("Tp", Tp, units="N/m")
ivc.add_output("Np", Np, units="N/m")
ivc.add_output("chord", chord, units="m")
ivc.add_output("theta", theta, units="rad")
ivc.add_output("A", A, units="m**2")
ivc.add_output("Iyy", Iyy, units="m**4")
ivc.add_output("Izz", Izz, units="m**4")
ivc.add_output("Iyz", Iyz, units="m**4")
prob.model.add_subsystem("ivc", ivc, promotes=["*"])

solver_comp = make_component(SolverComp(rho=2780.0, E=72.4e9, nu=0.33, Rhub=Rhub, span=span, nelems=nelems))
prob.model.add_subsystem("solver_comp", solver_comp, promotes=["*"])

# prob.setup()
# prob.run_model()
# prob.check_partials(compact_print=true)

include("mass_component.jl")

prob = om.Problem()
ivc = om.IndepVarComp()
ivc.add_output("A", A, units="m**2")

solver_comp = make_component(MassComp(rho=2600.0, span=span, nelems=nelems))
prob.model.add_subsystem("mass_comp", solver_comp, promotes=["*"])

# prob.setup()
# prob.run_model()
# prob.check_partials(compact_print=true)

include("von_mises_component.jl")

num_stress_eval_points = 20
sigma1 = collect(range(-1.5, 1.5, length=2*nelems*num_stress_eval_points))

prob = om.Problem()
ivc = om.IndepVarComp()
ivc.add_output("sigma1", sigma1)

vm_comp = make_component(VonMisesComp(nelems=nelems, num_stress_eval_points=num_stress_eval_points))
prob.model.add_subsystem("vm_comp", vm_comp, promotes=["*"])

prob.setup()
prob.run_model()
prob.check_partials(compact_print=true)