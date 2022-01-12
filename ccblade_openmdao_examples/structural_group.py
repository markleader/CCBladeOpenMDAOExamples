import numpy as np

import openmdao.api as om
from julia.OpenMDAO import make_component
import julia.Main as julia
from ccblade_openmdao_examples.gxbeam_openmdao_component import AreaComp
from ccblade_openmdao_examples.gxbeam_openmdao_component import SolverComp
from ccblade_openmdao_examples.gxbeam_openmdao_component import StressComp
from ccblade_openmdao_examples.gxbeam_openmdao_component import MassComp

class StructuralGroup(om.Group):
    def initialize(self):
        self.options.declare("nelems", types=int)
        self.options.declare("num_stress_eval_points", types=int)
        self.options.declare("span", default=12.0*0.0254)
        self.options.declare("Rhub", default=2.4*0.0254)
        self.options.declare("ys", default=345e6)

        return

    def setup(self):

        num_stress_eval_points = self.options["num_stress_eval_points"]
        nelems = self.options["nelems"]
        span = self.options["span"]
        Rhub = self.options["Rhub"]
        ys = self.options["ys"]

        area_comp = make_component(
            AreaComp(
                nelems=nelems, A_ref=821.8, Iyy_ref=23543.4, Izz_ref=5100.8))

        solver_comp = make_component(
            SolverComp(
                rho=2780.0, E=72.4e9, nu=0.33, Rhub=Rhub, span=span, nelems=nelems))

        stress_comp = make_component(
            StressComp(
                nelems=nelems, num_stress_eval_points=num_stress_eval_points, ys=ys))

        mass_comp = make_component(
            MassComp(rho=2780.0, span=span, nelems=nelems))

        self.add_subsystem(name="area_comp",
                           subsys=area_comp,
                           promotes_inputs=["*"],
                           promotes_outputs=["*"])

        self.add_subsystem(name="solver_comp",
                           subsys=solver_comp,
                           promotes_inputs=["*"])

        self.add_subsystem(name="stress_comp",
                           subsys=stress_comp,
                           promotes_inputs=["chord", "theta", "A", "Iyy", "Izz", "Iyz"],
                           promotes_outputs=["sigma1"])

        self.add_subsystem(name="mass_comp",
                           subsys=mass_comp,
                           promotes_inputs=["A"],
                           promotes_outputs=["m"])

        self.connect("solver_comp.Fx", "stress_comp.Fx")
        self.connect("solver_comp.My", "stress_comp.My")
        self.connect("solver_comp.Mz", "stress_comp.Mz")

        return