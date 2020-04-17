
import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])


import CompuCellSetup


sim,simthread = CompuCellSetup.getCoreSimulationObjects()
            
CompuCellSetup.initializeSimulationObjects(sim,simthread)
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()
        
from BuildWall3DSteppables import ConstraintInitializerSteppableAdder
steppableInstance=ConstraintInitializerSteppableAdder(sim,_frequency=1)
steppableRegistry.registerSteppable(steppableInstance)

from BuildWall3DSteppables import GrowthSteppableLinear
steppableInstance=GrowthSteppableLinear(sim,_frequency=10)
steppableRegistry.registerSteppable(steppableInstance)

from BuildWall3DSteppables import MitosisSteppableAdder
steppableInstance=MitosisSteppableAdder(sim,_frequency=10)
steppableRegistry.registerSteppable(steppableInstance)

from BuildWall3DSteppables import BuildWall3DSteppable
steppableInstance=BuildWall3DSteppable(sim,_frequency=10)
steppableRegistry.registerSteppable(steppableInstance)
        
from BuildWall3DSteppables import DeathSteppable
steppableInstance=DeathSteppable(sim,_frequency=5)
steppableRegistry.registerSteppable(steppableInstance)       

from BuildWall3DSteppables import neighbourdata
steppableInstance=neighbourdata(sim,_frequency=10)
steppableRegistry.registerSteppable(steppableInstance) 

from BuildWall3DSteppables import tracking
steppableInstance=tracking(sim,_frequency=10)
steppableRegistry.registerSteppable(steppableInstance) 
        
CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        