import numpy as np

import openmdao.api as om
from julia.OpenMDAO import make_component
import julia.Main as julia
from ccblade_openmdao_examples.gxbeam_openmdao_component import AreaComp
from ccblade_openmdao_examples.gxbeam_openmdao_component import SolverComp
from ccblade_openmdao_examples.gxbeam_openmdao_component import StressComp
from ccblade_openmdao_examples.gxbeam_openmdao_component import VonMisesComp
from ccblade_openmdao_examples.gxbeam_openmdao_component import MassComp

class StructuralGroup(om.Group):
    def initialize(self):
        self.options.declare("nelems", types=int, desc="number of elements - same for aerodynamic and structural analyses")
        self.options.declare("num_stress_eval_points", types=int, desc="number of points to query stress per element")
        self.options.declare("span", default=12.0*0.0254, desc="distance to tip of the blade (m)")
        self.options.declare("Rhub", default=2.4*0.0254, desc="hub radius (m)")
        self.options.declare("ys", default=345e6, desc="material yield stress (Pa)")
        self.options.declare("rho", default=2780.0, desc="material density (kg/m^3)")
        self.options.declare("zrel_cm", default=0.41, lower=0.0, upper=1.0, desc="distance from the airfoil tip to the center of mass, as a percent of chord")
        self.options.declare("zrel_rot", default=0.25, lower=0.0, upper=1.0, desc="distance from the airfoil tip to the axis of twist, as a percent of chord")

        return

    def setup(self):

        num_stress_eval_points = self.options["num_stress_eval_points"]
        nelems = self.options["nelems"]
        span = self.options["span"]
        Rhub = self.options["Rhub"]
        ys = self.options["ys"]
        rho = self.options["rho"]
        zrel_cm = self.options["zrel_cm"]
        zrel_rot = self.options["zrel_rot"]

        area_comp = make_component(
            AreaComp(
                nelems=nelems, A_ref=821.8, Iyy_ref=23543.4, Izz_ref=5100.8, zrel_cm=zrel_cm, zrel_rot=zrel_rot))

        solver_comp = make_component(
            SolverComp(
                rho=2780.0, E=72.4e9, nu=0.33, Rhub=Rhub, span=span, nelems=nelems))

        stress_comp = make_component(
            StressComp(
                nelems=nelems, num_stress_eval_points=num_stress_eval_points, ys=ys, zrel_rot=zrel_rot))

        vm_comp = make_component(
            VonMisesComp(
                nelems=nelems, num_stress_eval_points=num_stress_eval_points, k0=0.01))

        mass_comp = make_component(
            MassComp(rho=rho, span=span, nelems=nelems))

        self.add_subsystem(name="area_comp",
                           subsys=area_comp,
                           promotes_inputs=["*"],
                           promotes_outputs=["*"])

        self.add_subsystem(name="solver_comp",
                           subsys=solver_comp,
                           promotes_inputs=["*"])

        self.add_subsystem(name="stress_comp",
                           subsys=stress_comp,
                           promotes_inputs=["chord", "theta", "A", "Iyy", "Izz", "Iyz"])

        self.add_subsystem(name="vm_comp",
                           subsys=vm_comp,
                           promotes_outputs=["sigma_vm"])

        self.add_subsystem(name="mass_comp",
                           subsys=mass_comp,
                           promotes_inputs=["A"],
                           promotes_outputs=["m"])

        self.connect("solver_comp.Fx", "stress_comp.Fx")
        self.connect("solver_comp.My", "stress_comp.My")
        self.connect("solver_comp.Mz", "stress_comp.Mz")
        self.connect("stress_comp.sigma1", "vm_comp.sigma1")

        return