import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])


import CompuCellSetup


sim,simthread = CompuCellSetup.getCoreSimulationObjects()
            
# add extra attributes here
            
CompuCellSetup.initializeSimulationObjects(sim,simthread)
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()
        

from CompetitionSteppables import ConstraintInitializerSteppableAdder
ConstraintInitializerSteppableAdderInstance=ConstraintInitializerSteppableAdder(sim,_frequency=1)
steppableRegistry.registerSteppable(ConstraintInitializerSteppableAdderInstance)

from CompetitionSteppables import GrowthSteppableLinear
GrowthSteppableLinearInstance=GrowthSteppableLinear(sim,_frequency=10)
steppableRegistry.registerSteppable(GrowthSteppableLinearInstance)


from CompetitionSteppables import MitosisSteppableAdder
MitosisSteppableAdderInstance=MitosisSteppableAdder(sim,_frequency=10)
steppableRegistry.registerSteppable(MitosisSteppableAdderInstance)

from CompetitionSteppables import CellMotilitySteppable
CellMotilitySteppableInstance=CellMotilitySteppable(sim,_frequency=10)
steppableRegistry.registerSteppable(CellMotilitySteppableInstance)

from CompetitionSteppables import DeathSteppable
DeathSteppableInstance=DeathSteppable(sim,_frequency=10)
steppableRegistry.registerSteppable(DeathSteppableInstance)

# from CompetitionSteppables import DeathSteppable_Elasticity
# DeathSteppable_ElasticityInstance=DeathSteppable_Elasticity(sim,_frequency=10)
# steppableRegistry.registerSteppable(DeathSteppable_ElasticityInstance)

# from CompetitionSteppables import DeathSteppablePerimiter
# DeathSteppablePerimiterInstance=DeathSteppablePerimiter(sim,_frequency=10)
# steppableRegistry.registerSteppable(DeathSteppablePerimiterInstance)

from CompetitionSteppables import neighbourdata
neighbourdataInstance=neighbourdata(sim,_frequency=10)
steppableRegistry.registerSteppable(neighbourdataInstance) 

from CompetitionSteppables import tracking
trackingInstance=tracking(sim,_frequency=10)
steppableRegistry.registerSteppable(trackingInstance)
 
from CompetitionSteppables import cleanup
cleanupInstance=cleanup(sim,_frequency=1000)
steppableRegistry.registerSteppable(cleanupInstance)





        
CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        