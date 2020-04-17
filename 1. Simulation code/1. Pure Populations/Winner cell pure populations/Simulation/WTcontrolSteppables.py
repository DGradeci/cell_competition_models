# Steppables & Cellular Rules 
"""
By Daniel Gradeci
London Center of Nanotechnolgy
UCL Physics & Astronomy

"""

from __future__ import division
from PlayerPython import * 
import CompuCellSetup
from PySteppables import *
import CompuCell
import sys
import numpy as np
import random
from datetime import datetime
from tempfile import TemporaryFile


from PySteppablesExamples import MitosisSteppableBase
            
            
global minvol
global maxvol
global relaxtime
global adderlist
adderlist=[]
growthratewt=6.0# growth rate 
CI=0.1 # Contact Inhibition Sensitvity 
minvol=1200
maxvol=2200
relaxtime=200
          
            
            
#Initiate Target Volume, with some noise and Lambda
class ConstraintInitializerSteppable(SteppableBasePy):
    
#Initiate Target Volume, with some noise and Lambda 
#for an adder growth rule
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        for cell in self.cellList:
            if cell.type==1 or cell.type==2 or cell.type==3 or cell.type==4:
                cell.targetVolume=random.normalvariate(1800, 500)
                round(cell.targetVolume)        
                cell.lambdaVolume=2.0
           
        
        


class ConstraintInitializerSteppableAdder(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        deletecells=random.randint(1,50)
        deletecells1=list(range(deletecells))
        for cell in self.cellList:
            if cell.id in deletecells1:
                cell.type=5
            x_shift=random.randint(0, 1200)-int(round(cell.xCOM));
            y_shift=random.randint(0, 1600)-int(round(cell.yCOM));
            self.moveCell(cell,[x_shift,y_shift,0])
            if cell.type==1 or cell.type==2 or cell.type==3 or cell.type==4:
                cell.targetVolume=random.uniform(800, 3000)
                round(cell.targetVolume)          
                cell.lambdaVolume=2.0
                id=cell.id   
                add1=[id, 1800,cell.type]
                adderlist.append(add1)
#                 volList.append(cell.targetVolume)  
   






#Growth of cells for mitosis, with noise, also change lambda with time so size reduced.

#cell growth with linear growth
class GrowthSteppableLinear(SteppableBasePy):
    def __init__(self,_simulator,_frequency=5):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        if mcs>relaxtime:
            for cell in self.cellList:
                if cell.type!=5:
                    cell.targetVolume+=round(random.normalvariate(growthratewt, 2.5))*np.exp(-(CI/mcs)*((cell.volume-cell.targetVolume)**2))  # a = softness parameter
#                     tvadded+=round(random.normalvariate(growthratewt, 2.5))*np.exp(-(CI/mcs)*((cell.volume-cell.targetVolume)**2))       
                else:
                    cell.targetVolume=0    
    # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        
        # field=CompuCell.getConcentrationField(self.simulator,"PUT_NAME_OF_CHEMICAL_FIELD_HERE")
        # pt=CompuCell.Point3D()
        # for cell in self.cellList:
            # pt.x=int(cell.xCOM)
            # pt.y=int(cell.yCOM)
            # pt.z=int(cell.zCOM)
            # concentrationAtCOM=field.get(pt)
            # cell.targetVolume+=0.01*concentrationAtCOM  # you can use here any fcn of concentrationAtCOM     
        
        









#cell growth with exponential growth
class GrowthSteppableExp(SteppableBasePy):
    
#division for adder mechanism    
    
    
    
    
    def __init__(self,_simulator,_frequency=5):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        if mcs>relaxtime:
            for cell in self.cellList:
                if cell.type!=5:
#                     cell.targetVolume+=cell.targetVolume*np.exp(0.00000000035) 
                      cell.targetVolume=cell.targetVolume*np.exp(random.uniform(0.0015,0.0020)) 
                else:
                    cell.targetVolume=0    
    # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        
        # field=CompuCell.getConcentrationField(self.simulator,"PUT_NAME_OF_CHEMICAL_FIELD_HERE")
        # pt=CompuCell.Point3D()
        # for cell in self.cellList:
            # pt.x=int(cell.xCOM)
            # pt.y=int(cell.yCOM)
            # pt.z=int(cell.zCOM)
            # concentrationAtCOM=field.get(pt)
            # cell.targetVolume+=0.01*concentrationAtCOM  # you can use here any fcn of concentrationAtCOM     
        
        


class MitosisSteppableAdder(MitosisSteppableBase):
    
#division for sizer mechanism    
    
    
    def __init__(self,_simulator,_frequency=10):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
        
    def start(self):
        global parentlist 
        global parentid  
        global childid 
        global Volumeanalysis
#         global adderlist
        parentlist=[]
        parentid=0
        childid=0
        Volumeanalysis=[]
#         adderlist=[]
#         for cell in self.cellList:
#             id=cell.id
#             birthvol=cell.targetVolume
#             add1=[id,birthvol]
#             adderlist.append(add1)
    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        if mcs>relaxtime :
            for cell in self.cellList:
                for x in adderlist:
                    if x[0]==cell.id and cell.volume-x[1] > (1800-(mcs/60)) and (cell.type==1 or cell.type==2 or cell.type==3 or cell.type==4 ):
                        cells_to_divide.append(cell)  
                        Volumeanalysis.append([x[0],(mcs-relaxtime)/float(10),x[1],cell.volume,cell.type])                       
        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
            # self.divideCellRandomOrientation(cell)
            # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
            # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
            self.divideCellAlongMinorAxis(cell)                               # this is a valid option
 
            parentid=self.parentCell.id
            childid=self.childCell.id
            par1=[parentid,childid,(mcs-relaxtime)/float(10)]
            parentlist.append(par1)
            birthvolchild=self.childCell.targetVolume
            addchild=[childid,birthvolchild]
            adderlist.append(addchild)
#             for x in adderlist:
#                 if x[0]==cell.id:
#                     x[1]/=2


#     def storeparents(self,mcs):  
#         for cell in self.cellList:
#             if cell.type !=1:   
#                 parentid=self.parentCell.id
#                 par1=[parentid,mcs]
#                 parentlist.append(par1)
    def updateAttributes(self):
        self.parentCell.targetVolume /=2 # reducing parent target volume             
        self.cloneParent2Child() 

 

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.cloneAttributes(sourceCell=self.parentCell, targetCell=self.childCell, no_clone_key_dict_list = [attrib1, attrib2] )
        
        
        if self.parentCell.type==1:
            self.childCell.type=2
            self.parentCell.type=2
        elif self.parentCell.type==2:
            self.childCell.type=3
            self.parentCell.type=3
        elif self.parentCell.type==3:
            self.childCell.type=4   
            self.parentCell.type=4
        elif self.parentCell.type==4:
            self.childCell.type=1   
            self.parentCell.type=1
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\\parentlist-%s'%datetime.now().strftime('%H-%M-%m-%d'), parentlist)   
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\\Volumeanalysis-%s'%datetime.now().strftime('%H-%M-%m-%d'),Volumeanalysis)  
        
        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.cloneAttributes(sourceCell=self.parentCell, targetCell=self.childCell, no_clone_key_dict_list = [attrib1, attrib2] )
        
        
     
        


class MitosisSteppable(MitosisSteppableBase):
    
#division for timer mechanism 
    def __init__(self,_simulator,_frequency=10):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
        
    def start(self):
        global parentlist 
        global parentid  
        global childid 
        parentlist=[]
        parentid=0
        childid=0
    
    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        if mcs>relaxtime :
            for cell in self.cellList:
                if cell.id%2==random.randint(0,1):
                    if cell.type==1 or cell.type==2 or cell.type==3 :
                        if cell.volume>random.randint(minvol*2,maxvol*1.5):
                            cells_to_divide.append(cell)
                    
        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
            self.divideCellRandomOrientation(cell)
            # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
            # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
            # self.divideCellAlongMinorAxis(cell)                               # this is a valid option
 
            parentid=self.parentCell.id
            childid=self.childCell.id
            
            par1=[parentid,childid,mcs]
            parentlist.append(par1)


#     def storeparents(self,mcs):  
#         for cell in self.cellList:
#             if cell.type !=1:   
#                 parentid=self.parentCell.id
#                 par1=[parentid,mcs]
#                 parentlist.append(par1)
    def updateAttributes(self):
        self.parentCell.targetVolume /= 2 # reducing parent target volume             
        self.cloneParent2Child() 

 
        
        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.cloneAttributes(sourceCell=self.parentCell, targetCell=self.childCell, no_clone_key_dict_list = [attrib1, attrib2] )
        
        
        if self.parentCell.type==1:
            self.childCell.type=2
            self.parentCell.type=2
        elif self.parentCell.type==2:
            self.childCell.type=3
            self.parentCell.type=3
        elif self.parentCell.type==3:
            self.childCell.type=4   
            self.parentCell.type=4
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\parentlist-%s'%datetime.now().strftime('%H-%M-%m-%d'), parentlist)         
        
        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.cloneAttributes(sourceCell=self.parentCell, targetCell=self.childCell, no_clone_key_dict_list = [attrib1, attrib2] )
        
        
     
        








class MitosisSteppableTimer(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=10):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
        
    def start(self):
        global parentlist 
        global parentid  
        global childid 
        global adderlist
        parentlist=[]
        parentid=0
        childid=0
        adderlist=[]
        for cell in self.cellList:
            id=cell.id
            celltimer=random.randint(-1500,0)
            timer1=[id,celltimer]
            adderlist.append(timer1)
    def step(self,mcs):
        # print "INSIDE MITOSIS STEPPABLE"
        cells_to_divide=[]
        if mcs>relaxtime :
            for cell in self.cellList:
                for x in adderlist:
                    if x[0]==cell.id and mcs-x[1] > 2500 and (cell.type==1 or cell.type==2 or cell.type==3 or cell.type==4):
                        cells_to_divide.append(cell) 
                        x[1]=mcs                 
        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
            self.divideCellRandomOrientation(cell)
            # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
            # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
            # self.divideCellAlongMinorAxis(cell)                               # this is a valid option
 
            parentid=self.parentCell.id
            childid=self.childCell.id
            par1=[parentid,childid,mcs]
            parentlist.append(par1)
            birthvolchild=self.childCell.targetVolume
            addchild=[childid,mcs]
            adderlist.append(addchild)


#     def storeparents(self,mcs):  
#         for cell in self.cellList:
#             if cell.type !=1:   
#                 parentid=self.parentCell.id
#                 par1=[parentid,mcs]
#                 parentlist.append(par1)
    def updateAttributes(self):
        self.parentCell.targetVolume /=2 # reducing parent target volume             
        self.cloneParent2Child() 

 

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.cloneAttributes(sourceCell=self.parentCell, targetCell=self.childCell, no_clone_key_dict_list = [attrib1, attrib2] )
        
        
        if self.parentCell.type==1:
            self.childCell.type=2
            self.parentCell.type=2
        elif self.parentCell.type==2:
            self.childCell.type=3
            self.parentCell.type=3
        elif self.parentCell.type==3:
            self.childCell.type=4   
            self.parentCell.type=4
        elif self.parentCell.type==4:
            self.childCell.type=1   
            self.parentCell.type=1
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\\parentlist-%s'%datetime.now().strftime('%H-%M-%m-%d'), parentlist)         
        
        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.cloneAttributes(sourceCell=self.parentCell, targetCell=self.childCell, no_clone_key_dict_list = [attrib1, attrib2] )
        
        
     
        










#Density Dependent death.
class DeathSteppable(SteppableBasePy):
    
#plotting cell counts v time.    
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
# #        self.pW=self.addNewPlotWindow(_title='Appoptosis Count',_xAxisTitle='Time (Frames)',_yAxisTitle='Apoptosis Count',_xScaleType='linear',_yScaleType='log')
#         self.pW=self.addNewPlotWindow(_title='Density Count',_xAxisTitle='Density',_yAxisTitle='Apoptosis Count',_xScaleType='linear',_yScaleType='log')
# #        self.pW.addPlot('Apoptosis Count',_style='Lines',_color='black',_size=5)  
#         self.pW.addPlot('Density Count',_style='Dots',_color='green',_size=5)
        self.pW=self.addNewPlotWindow(_title='Local Density v Time',_xAxisTitle='Time(Frames)',_yAxisTitle='Avg Local Density',_xScaleType='linear',_yScaleType='linear')
        self.pW.addPlot('ld WT',_style='Dots',_color='red',_size=5)  
        self.pW.addPlot('p_apo_wt',_style='Dots',_color='blue',_size=5)    

    def step(self,mcs):
        deathcount=0
        dens=0
        totdens=0
        cellcount=0
        
        if mcs>relaxtime:
            for cell in self.cellList:   
                if cell.type!=5:
                    if (cell.xCOM>200 or cell.xCOM<1000) and (cell.yCOM>200 or cell.xCOM<1400):
                        cellneighbours=[]
                        nieghbourvolume=[]
                        numberOfCellNeighbours=[]
                        dens=(1/cell.volume)
                        for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                            if neighbor:
                                cellneighbours.append(neighbor.id)   
                                nieghbourvolume.append(neighbor.volume)
                                dens+=(1/neighbor.volume)
                        cn= (str(w) for w in cellneighbours)
                        den= (str(k) for k in nieghbourvolume)
                        li=list(nieghbourvolume)
                        ne=len(cellneighbours)
                        totdens+=dens
                        cellcount+=1
                        if cell.type==1 or cell.type==2 or cell.type==3 or cell.type==4:
                            if dens >0 and 0.00002>random.random():
                                cell.type=5
                                cell.targetVolume=0
                                cell.lambdaVolume=2
                                deathcount+=1
            time = (mcs-relaxtime)/float(10)              
            self.pW.addDataPoint('ld WT',time,totdens/cellcount)  
            self.pW.savePlotAsData('Local_Denstiy_evolution_wt.txt')   
            self.pW.addDataPoint('p_apo_wt',totdens/cellcount,deathcount/cellcount)  
            self.pW.savePlotAsData('p_apo_wt.txt')             
            
#             self.pW.addDataPoint('Apoptosis Count',(mcs/float(10)),deathcount) # name of the data series, x, y    
#             self.pW.addDataPoint('Density Count',dens,deathcount)   
    def finish(self):
        return
        
        
class neighbourdata(SteppableBasePy):
    
#storing information on cells.       
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)   
        
    def start(self):
        self.pW=self.addNewPlotWindow(
        _title='CellCount v Time',
        _xAxisTitle='Time(Frames)',
        _yAxisTitle='Cell Count',
        _xScaleType='linear',
        _yScaleType='linear')
        
        self.pW.addPlot('cellcount WT',_style='Dots',_color='green',_size=5)   
#         self.pW=self.addNewPlotWindow(_title='Cellvolume v Time',_xAxisTitle='Time(Frames)',_yAxisTitle='Cell volume',_xScaleType='linear',_yScaleType='log')
#         self.pW.addPlot('cellvol WT',_style='Dots',_color='blue',_size=5) 
    def step(self,mcs):
        cellcount=0
        cellvolume=0
        for cell in self.cellList:
            if cell.type != 5:
                if cell.type==1 or cell.type==2 or cell.type==3 or cell.type==4:
                    cellvolume+=cell.volume 
                    cellcount+=1
        time = (mcs-relaxtime)/float(10)       
        self.pW.addDataPoint('cellcount WT',time,cellcount)  
        self.pW.savePlotAsData('Cell_Count.txt')
#         self.pW.addDataPoint("cellvol  WT",time,cellvolume/cellcount) 
    def finish(self):
        return
 




class tracking(SteppableBasePy): 
  
# clean up dead/extruded cells  
    
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)     
    def start(self): 
        self.pW=self.addNewPlotWindow(_title='CellVol v Time',_xAxisTitle='Time(Frames)',_yAxisTitle='Cell Count',_xScaleType='linear',_yScaleType='log')
        self.pW.addPlot('cellvol',_style='Dots',_color='red',_size=5) 
        global trackingfile
        trackingfile=[]
    def step(self,mcs): 
        cellcount=0
        cellvolume=0
        for cell in self.cellList:
            if cell.type==5 :
                state=1
            else:
                state=0
            if cell.type==1 or cell.type==2 or cell.type==3 or cell.type==4:
                cellvolume+=cell.volume 
                cellcount+=1
            time = (mcs-relaxtime)/float(10)      
            ar1= [cell.xCOM,cell.yCOM,int(time),int(cell.id),int(cell.type),state]
            trackingfile.append(ar1)
        self.pW.addDataPoint("cellvol",time,cellvolume/cellcount)   
        self.pW.savePlotAsData('Cell_Vol.txt')  
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\\tracking-%s'%datetime.now().strftime('%H-%M-%m-%d'), trackingfile)
        return


 
class cleanup(SteppableBasePy):   
    

#stored Local density information    
        
    def __init__(self,_simulator,_frequency=100):
        SteppableBasePy.__init__(self,_simulator,_frequency)     
    def start(self): 
        self.pW=self.addNewPlotWindow(_title='extrusions',_xAxisTitle='time',_yAxisTitle='extrusion',_xScaleType='linear',_yScaleType='linear')
        self.pW.addPlot('ext',_style='Dots',_color='red',_size=5)
    def step(self,mcs): 
        if mcs>relaxtime:
            ext=0
            for cell in self.cellList:
                if cell.volume<500:
                    ext+=1
                    cell.type=5
            time = (mcs-relaxtime)/float(10)          
            self.pW.addDataPoint("ext",time,ext)   
            self.pW.savePlotAsData('extrusions.txt')  
                    
                












class P_locdens(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        self.pW=self.addNewPlotWindow(_title='P_div',_xAxisTitle='Avg Local Density',_yAxisTitle='P_div',_xScaleType='linear',_yScaleType='linear')
        self.pW.addPlot('p_div_wt',_style='Dots',_color='blue',_size=5)   

    def step(self,mcs):
        p_div=0
        totdens=0
        cellcount=0
        dens=0
        Div_list=len(parentlist)
        if mcs>relaxtime:    
            for cell in self.cellList:   
                if (cell.xCOM>200 or cell.xCOM<1000) and (cell.yCOM>200 or cell.xCOM<1400):
                    cellneighbours=[]
                    nieghbourvolume=[]
                    numberOfCellNeighbours=[]
                    dens=(1/cell.volume)
                    for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                        if neighbor:
                            cellneighbours.append(neighbor.id)   
                            nieghbourvolume.append(neighbor.volume)
                            dens+=(1/neighbor.volume)
                    li=list(nieghbourvolume)
                    totdens+=dens
                    cellcount+=1             
#             self.pW.addDataPoint('p_div_wt',totdens/cellcount,Div_list/cellcount)  
            self.pW.savePlotAsData('p_div_wt.txt')   