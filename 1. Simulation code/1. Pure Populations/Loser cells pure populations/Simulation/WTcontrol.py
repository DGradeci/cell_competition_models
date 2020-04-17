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
        

# from WTcontrolSteppables import ConstraintInitializerSteppable
# ConstraintInitializerSteppableInstance=ConstraintInitializerSteppable(sim,_frequency=1)
# steppableRegistry.registerSteppable(ConstraintInitializerSteppableInstance)
        


from WTcontrolSteppables import ConstraintInitializerSteppableAdder
ConstraintInitializerSteppableAdderInstance=ConstraintInitializerSteppableAdder(sim,_frequency=1)
steppableRegistry.registerSteppable(ConstraintInitializerSteppableAdderInstance)
        


from WTcontrolSteppables import GrowthSteppableLinear
GrowthSteppableLinearInstance=GrowthSteppableLinear(sim,_frequency=10)
steppableRegistry.registerSteppable(GrowthSteppableLinearInstance)

# from WTcontrolSteppables import GrowthSteppableExp
# GrowthSteppableExpInstance=GrowthSteppableExp(sim,_frequency=5)
# steppableRegistry.registerSteppable(GrowthSteppableExpInstance)        

# from WTcontrolSteppables import MitosisSteppable
# MitosisSteppableInstance=MitosisSteppable(sim,_frequency=20)
# steppableRegistry.registerSteppable(MitosisSteppableInstance)

from WTcontrolSteppables import MitosisSteppableAdder
MitosisSteppableAdderInstance=MitosisSteppableAdder(sim,_frequency=10)
steppableRegistry.registerSteppable(MitosisSteppableAdderInstance)

# from WTcontrolSteppables import MitosisSteppableTimer
# MitosisSteppableTimerInstance=MitosisSteppableTimer(sim,_frequency=10)
# steppableRegistry.registerSteppable(MitosisSteppableTimerInstance)
        


from WTcontrolSteppables import CellMotilitySteppable
CellMotilitySteppableInstance=CellMotilitySteppable(sim,_frequency=10)
steppableRegistry.registerSteppable(CellMotilitySteppableInstance)



from WTcontrolSteppables import DeathSteppable
DeathSteppableInstance=DeathSteppable(sim,_frequency=10)
steppableRegistry.registerSteppable(DeathSteppableInstance)

from WTcontrolSteppables import neighbourdata
neighbourdataInstance=neighbourdata(sim,_frequency=10)
steppableRegistry.registerSteppable(neighbourdataInstance) 

from WTcontrolSteppables import tracking
trackingInstance=tracking(sim,_frequency=10)
steppableRegistry.registerSteppable(trackingInstance)
 
from WTcontrolSteppables import cleanup
cleanupInstance=cleanup(sim,_frequency=100)
steppableRegistry.registerSteppable(cleanupInstance)

from WTcontrolSteppables import P_locdens
P_locdensInstance=P_locdens(sim,_frequency=10)
steppableRegistry.registerSteppable(P_locdensInstance)


# from WTcontrolSteppables import ContactInhibition
# ContactInhibitionInstance=ContactInhibition(sim,_frequency=10)
# steppableRegistry.registerSteppable(ContactInhibitionInstance)

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        