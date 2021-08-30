using OpenMDAO: om, make_component

include("area_component.jl")

prob = om.Problem()

num_elem = 5
chord = collect(range(1.5*0.0254, 0.5*0.0254, length=num_elem))
twist = collect(range(1.0, 0.0, length=num_elem))

A_ref = 821.8
Iyy_ref = 23543.4
Izz_ref = 5100.8

ivc = om.IndepVarComp()
ivc.add_output("chord", chord, units="m")
ivc.add_output("twist", twist, units="rad")
prob.model.add_subsystem("ivc", ivc, promotes=["*"])

area_comp = make_component(AreaComp(num_elem=num_elem, A_ref=A_ref, Iyy_ref=Iyy_ref, Izz_ref=Izz_ref))  # Need to convert Julia obj to Python obj
prob.model.add_subsystem("area_comp", area_comp, promotes=["*"])

prob.setup()
prob.run_model()
#prob.check_partials(compact_print=true)

include("stress_component.jl")

prob = om.Problem()

num_elem = 5
Nx = collect(range(1000.0, 0.0, length=num_elem))
My = collect(range(100.0, 0.0, length=num_elem))
Mz = collect(range(500.0, 0.0, length=num_elem))
chord = collect(range(1.5*0.0254, 0.5*0.0254, length=num_elem))
twist = collect(range(1.0, 0.0, length=num_elem))
A = collect(range(0.01, 0.005, length=num_elem))
Iyy = collect(range(0.002, 0.001, length=num_elem))
Izz = collect(range(0.002, 0.001, length=num_elem))
Iyz = collect(range(0.001, 0.0, length=num_elem))

ivc = om.IndepVarComp()
ivc.add_output("Nx", Nx, units="N")
ivc.add_output("My", My, units="N*m")
ivc.add_output("Mz", Mz, units="N*m")
ivc.add_output("chord", chord, units="m")
ivc.add_output("twist", twist, units="rad")
ivc.add_output("A", A, units="m**2")
ivc.add_output("Iyy", Iyy, units="m**4")
ivc.add_output("Izz", Izz, units="m**4")
ivc.add_output("Iyz", Iyz, units="m**4")
prob.model.add_subsystem("ivc", ivc, promotes=["*"])

area_comp = make_component(StressComp(num_elem=num_elem, num_stress_eval_points=20))
prob.model.add_subsystem("stress_comp", area_comp, promotes=["*"])

prob.setup()
prob.run_model()
#prob.check_partials(compact_print=true)

include("gxbeam_component.jl")

prob = om.Problem()

span = 0.3
num_elem = 10
omega = 7110.0*2*pi/60
x = collect(range(0.06, span, length=num_elem+1))
Tp = collect(range(10.0, 50.0, length=num_elem+1))
Np = collect(range(10.0, 200.0, length=num_elem+1))
chord = collect(range(1.5*0.0254, 0.5*0.0254, length=num_elem))
twist = collect(range(1.0, 0.0, length=num_elem))

s = sin.(twist)
c = cos.(twist)
s2 = s.^2
c2 = c.^2

k = chord/100  # scale factor
A = (k.^2)*821.8
Iyy0 = (k.^4)*23543.4
Izz0 = (k.^4)*5100.8
Iyy = @. c2 *Iyy0 + s2 *Izz0
Izz = @. s2 *Iyy0 + c2 *Izz0
Iyz = @. (Iyy0 - Izz0)*s*c
# A = collect(range(0.01, 0.005, length=num_elem))
# Iyy = collect(range(0.002, 0.001, length=num_elem))
# Izz = collect(range(0.002, 0.001, length=num_elem))
# Iyz = collect(range(0.001, 0.0, length=num_elem))

ivc = om.IndepVarComp()
ivc.add_output("omega", omega, units="rad/s")
ivc.add_output("Tp", Tp, units="N/m")
ivc.add_output("Np", Np, units="N/m")
ivc.add_output("chord", chord, units="m")
ivc.add_output("twist", twist, units="rad")
ivc.add_output("A", A, units="m**2")
ivc.add_output("Iyy", Iyy, units="m**4")
ivc.add_output("Izz", Izz, units="m**4")
ivc.add_output("Iyz", Iyz, units="m**4")
prob.model.add_subsystem("ivc", ivc, promotes=["*"])

solver_comp = make_component(SolverComp(rho=2600.0, E=70e9, nu=0.33, x=x))
prob.model.add_subsystem("solver_comp", solver_comp, promotes=["*"])

prob.setup()
prob.run_model()
prob.check_partials(compact_print=true)