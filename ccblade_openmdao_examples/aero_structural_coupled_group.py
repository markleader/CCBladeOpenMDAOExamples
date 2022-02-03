import numpy as np

import openmdao.api as om
from julia.OpenMDAO import make_component
import julia.Main as julia

from ccblade_openmdao_examples.ccblade_openmdao_component import BEMTRotorCAComp
from .structural_group import StructuralGroup

class AeroStructuralGroup(om.Group):
    def initialize(self):
        self.options.declare("nelems", types=int, desc="number of elements - same for aerodynamic and structural analyses")
        self.options.declare("num_stress_eval_points", types=int, desc="number of points to query stress per element")
        self.options.declare("Rtip", default=12.0*0.0254, desc="distance to tip of the blade (m)")
        self.options.declare("Rhub", default=2.4*0.0254, desc="hub radius (m)")
        self.options.declare("ys", default=345e6, desc="material yield stress (Pa)")
        self.options.declare("rho", default=2780.0, desc="material density (kg/m^3)")
        self.options.declare("zrel_cm", default=0.41, lower=0.0, upper=1.0, desc="distance from the airfoil tip to the center of mass, as a percent of chord")
        self.options.declare("zrel_rot", default=0.25, lower=0.0, upper=1.0, desc="distance from the airfoil tip to the axis of twist, as a percent of chord")

        self.options.declare("af_fname", types=str, desc="filename of airfoil data")
        self.options.declare("c", types=float, desc="constant chord (m)")
        self.options.declare("Re_exp", default=0.6)
        self.options.declare("num_operating_points", default=1)
        self.options.declare("num_blades", types=int)
        self.options.declare("rho0", types=float, desc="air density (kg/m^3)")
        self.options.declare("mu", types=float)
        self.options.declare("speefofsound", types=float, desc="(m/s)")

        return

    def setup(self):

        num_stress_eval_points = self.options["num_stress_eval_points"]
        nelems = self.options["nelems"]
        Rtip = self.options["Rtip"]
        Rhub = self.options["Rhub"]
        ys = self.options["ys"]
        rho = self.options["rho"]
        zrel_cm = self.options["zrel_cm"]
        zrel_rot = self.options["zrel_rot"]

        af_fname = self.options["af_fname"]
        c = self.options["c"]
        Re_exp = self.options["Re_exp"]
        num_operating_points = self.options["num_operating_points"]
        num_blades = self.options["num_blades"]
        rho0 = self.options["rho0"]
        mu = self.options["mu"]
        speedofsound = self.options["speedofsound"]

        struc_group = StructuralGroup(
            nelems=nelems, num_stress_eval_points=num_stress_eval_points,
            span=Rtip, Rhub=Rhub, ys=ys, rho=rho, zrel_cm=zrel_cm, zrel_rot=zrel_rot)

        aero_comp = make_component(
            BEMTRotorCAComp(
                af_fname=af_fname, cr75=c/Rtip, Re_exp=Re_exp,
                num_operating_points=num_operating_points, num_blades=num_blades,
                num_radial=nelems, rho=rho0, mu=mu, speedofsound=speedofsound))

        self.add_subsystem(name="struc_group", subsys=struc_group,
                           promotes_inputs=["omega", "Tp", "Np", "chord", "theta"],
                           promotes_outputs=["sigma_vm", "m"])

        self.add_subsystem(name="bemt_rotor_comp", subsys=aero_comp,
                           promotes_inputs=["Rhub", "Rtip", "radii", "chord", "theta", "v", "omega", "pitch"],
                           promotes_outputs=["thrust", "torque", "efficiency", "Tp", "Np"])

        self.linear_solver = om.DirectSolver()