import openmdao.api as om
from ccblade_openmdao_examples.structural_group import StructuralGroup

prob = om.Problem()

ivc = om.IndepVarComp()
ivc.add_output("Tp", val=Tp, units="N/m")
ivc.add_output("Np", val=Np, units="N/m")
ivc.add_output("chord", val=chord, units="m")
ivc.add_output("twist", val=twist, units="rad")
prob.model.add_subsystem("inputs_comp", ivc, promotes_outputs=["*"])