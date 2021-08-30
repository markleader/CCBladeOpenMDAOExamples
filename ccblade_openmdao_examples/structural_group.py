import numpy as np

import openmdao.api as om
from julia.OpenMDAO import make_component
import julia.Main as julia
from ccblade_openmdao_examples.area_component import AreaComp
from ccblade_openmdao_examples.gxbeam_component import GXBeamComp
from ccblade_openmdao_examples.stress_component import StressComp

class StructuralGroup(om.Group):
    def initialize(self):
        self.options.declare('num_elem', types=int)
        self.options.declare('num_radial', types=int)
        self.options.declare('num_stress_eval_points', types=int)

        return

    def setup(self):

        num_elem = self.options['num_elem']
        num_radial = self.options['num_radial']
        num_stress_eval_points = self.options['num_stress_eval_points']

        area_comp = make_component(
            AreaComp(
                num_elem=num_elem, A_ref=821.8, Iyy_ref=23543.4, Izz_ref=5100.8))

        gxbeam_comp = make_component(
            GXBeamComp(
                rho=2600.0, E=68.0e9, nu=0.33, span=12.0*0.0254, num_elem=num_elem, num_radial=num_radial))

        stress_comp = make_component(
            StressComp(
                num_stress_eval_points=num_stress_eval_points))

        self.add_subsystem(name='area_comp',
                           subsys=area_comp,
                           promotes_inputs=['*'],
                           promotes_outputs=['*'])

        self.add_subsystem(name='gxbeam_comp',
                           subsys=gxbeam_comp,
                           promotes_inputs=['*'])

        self.add_subsystem(name='stress_comp',
                           subsys=stress_comp,
                           promotes_inputs=['chord', 'twist', 'A', 'Iyy', 'Izz', 'Iyz'].
                           promotes_outputs=['sigma1'])

        self.connect('gxbeam_comp.Nx', 'stress_comp.Nx')
        self.connect('gxbeam_comp.My', 'stress_comp.My')
        self.connect('gxbeam_comp.Mz', 'stress_comp.Mz')

        return