
from cc3d import CompuCellSetup
        

from multicellularv2Steppables import AdhesionMoleculesSteppables

CompuCellSetup.register_steppable(steppable=AdhesionMoleculesSteppables(frequency=1))



from multicellularv2Steppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




from multicellularv2Steppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))




from multicellularv2Steppables import DeathSteppable

CompuCellSetup.register_steppable(steppable=DeathSteppable(frequency=1))


CompuCellSetup.run()
