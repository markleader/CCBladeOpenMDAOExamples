import os

from julia import Main

# All this is stolen from SciML/diffeqpy. I think after doing this, a user
# should have access to BEMTRotorComp.
script_dir = os.path.dirname(os.path.realpath(__file__))
Main.include(os.path.join(script_dir, "ccblade_openmdao_component.jl"))
Main.include(CCBlade)
BEMTRotorComp = Main.BEMTRotorComp
BEMTRotorCAComp = Main.BEMTRotorCAComp
