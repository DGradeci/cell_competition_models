from __future__ import division
from PlayerPython import * 
import CompuCellSetup
from PySteppables import *
from datetime import datetime
import CompuCell
import sys
import numpy as np
import random
from tempfile import TemporaryFile
from PySteppablesExamples import MitosisSteppableBase
        
    

# WT ROSENBLATT EXPERIMENT    

global minvol
global maxvol
global relaxtime
global adderlist
growthrate=3.0
adderlist=[]
# minvol=1200
# maxvol=3000
minvol=1100
maxvol=2200
relaxtime=200
          
            
            
#Initiate Target Volume, with some noise and Lambda
class ConstraintInitializerSteppableAdder(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        self.changeNumberOfWorkNodes(2) 
        volList = []
        self.pW=self.addNewPlotWindow(_title='initial seeding volumes' ,_xAxisTitle='initial volume',_yAxisTitle='N', _xScaleType='linear',_yScaleType='linear')
        self.pW.addHistogramPlot(_plotName='initialVol',_color='green',_alpha=100)# _alpha is transparency 0 is transparent, 255 is opaque        
        for cell in self.cellList:
#                 cell.targetVolume=random.randint(minvol-100,maxvol+100)    
#                 cell.targetVolume=random.normalvariate(2100, 300)
                cell.targetVolume=random.normalvariate(1600, 300)
                round(cell.targetVolume)
                cell.lambdaVolume=2.0
                id=cell.id   
                add1=[id,((maxvol+minvol)/2)-(cell.targetVolume/2.5),cell.type]
                adderlist.append(add1)
                volList.append(cell.targetVolume)  
        self.pW.addHistogram(plot_name='initialVol', value_array=volList, number_of_bins=100)     
        self.pW.savePlotAsData('initialVol.txt')    


class GrowthSteppableLinear(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        self.pW=self.addNewPlotWindow(
        _title='AvgVtAdded',
        _xAxisTitle='Time(Frames)',
        _yAxisTitle='AvgVtAdded',
        _xScaleType='linear',
        _yScaleType='linear')
        
        self.pW.addPlot('Avg Vt Added',_style='Dots',_color='green',_size=5) 
        self.pW.addPlot('Avg V',_style='Dots',_color='red',_size=5)        
    def step(self,mcs):
        cellcount=0
        tvadded=0
        totV=0
        if mcs>relaxtime:
            for cell in self.cellList:
                cellcount+=1
                totV+=cell.volume
                if cell.type!=5:
                    cell.targetVolume+=round(random.normalvariate(growthrate, 1.0))*np.exp(-(0.001/mcs)*((cell.volume-cell.targetVolume)**2))  # a = softness parameter
                    tvadded+=round(random.normalvariate(growthrate, 1.0))*np.exp(-(0.001/mcs)*((cell.volume-cell.targetVolume)**2))
                else:
                    cell.targetVolume=0  
        
        if cellcount==0:
            avgV=0
            avgtvadded=0
        else:
            avgV=totV/cellcount      
            avgtvadded=tvadded/cellcount  
        time = (mcs-relaxtime)/float(10)       
        self.pW.addDataPoint('Avg Vt Added',time,avgtvadded)  
        self.pW.savePlotAsData('avg_vt_added.txt')
        self.pW.addDataPoint('Avg V',time,avgV)  
        self.pW.savePlotAsData('avg_v.txt')
    # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        
        # field=CompuCell.getConcentrationField(self.simulator,"PUT_NAME_OF_CHEMICAL_FIELD_HERE")
        # pt=CompuCell.Point3D()
        # for cell in self.cellList:
            # pt.x=int(cell.xCOM)
            # pt.y=int(cell.yCOM)
            # pt.z=int(cell.zCOM)
            # concentrationAtCOM=field.get(pt)
            # cell.targetVolume+=0.01*concentrationAtCOM  # you can use here any fcn of concentrationAtCOM     
        
           
#division
class MitosisSteppableAdder(MitosisSteppableBase):  
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
                    if x[0]==cell.id and cell.volume-x[1] > 1800:
                        cells_to_divide.append(cell)  
                        Volumeanalysis.append([x[0],(mcs-relaxtime)/float(10),x[1],cell.volume,cell.type])                       
        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
            self.divideCellRandomOrientation(cell)
            # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
            # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
            # self.divideCellAlongMinorAxis(cell)                               # this is a valid option
 
            parentid=self.parentCell.id
            childid=self.childCell.id
            par1=[parentid,childid,(mcs-relaxtime)/float(10),cell.type]
            parentlist.append(par1)
            birthvolchild=self.childCell.targetVolume
            addchild=[childid,birthvolchild]
            adderlist.append(addchild)


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
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\parentlist-%s'%datetime.now().strftime('%H-%M-%m-%d'), parentlist)   
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\Volumeanalysis-%s'%datetime.now().strftime('%H-%M-%m-%d'),Volumeanalysis)  
        
        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.cloneAttributes(sourceCell=self.parentCell, targetCell=self.childCell, no_clone_key_dict_list = [attrib1, attrib2] )
        
        
     
        








class BuildWall3DSteppable(SteppableBasePy):    

    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        
        self.buildWall(self.WALL)

    def step(self,mcs):        

        if mcs<4800 or mcs>5000:
            for cell in self.cellList:
                
                cell.lambdaVecX=10.1*random.uniform(-1.0,1.0) # force component pointing along X axis 
                cell.lambdaVecY=10.1*random.uniform(-1.0,1.0) # force component pointing along Y axis 
                # print 'cell.lambdaVecX=',cell.lambdaVecX,' cell.lambdaVecY=',cell.lambdaVecY
        else:
            for cell in self.cellList:
                
#                 cell.lambdaVecY= (cell.yCOM-200)*2000 # force component pointing along Y axis 
                cell.lambdaVecX= (cell.xCOM-100)*2000# force component pointing along X axis 
        
#         if mcs==6500:
#             for cell in self.cellList:
#                     cell.targetVolume=cell.volume/1.5
#                     cell.lambdaVolume=1000            
        
        if mcs==5000:
#             for cell in self.cellList:
#                 cell.targetVolume=cell.volume*1.6
#                 cell.lambdaVolume=2.0
            xlist=[]
            ylist=[]
            for cell in self.cellList:
                if cell.type!=4:
                    xlist.append(cell.xCOM)  
                    ylist.append(cell.yCOM)
                
                
            xlat=max(xlist)+50
            print xlat
            print max(xlist)
            print max(ylist)
#             ylat=max(ylist)+50
            self.destroyWall()
            self.resizeAndShiftLattice(_newSize=(xlat,1000,1))
            self.buildWall(self.WALL)
        
        
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        




class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=5):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
       global P_apo
       P_apo=[]
    def step(self,mcs):
        deathcount=0
        dens=0
        totdens=0
        cellcount=0     
        if mcs>5000:
            for cell in self.cellList: 
                if cell.type!=5:
                    dens=(1/cell.volume)
                    if 0.00033074/(1 + np.exp(-235.8*(dens*9- 0.0152)))>random.random():
#                    if 0.00072194/(1 + np.exp(-509.4*(dens*9- 0.0026)))>random.random():#density dependent death 
                        cell.type=5
                        cell.targetVolume=0
                        cell.lambdaVolume=2
                        deathcount+=1
                    if cell.volume< 500:   #density dependent death  
                        cell.type=6
                        cell.targetVolume=0
                        cell.lambdaVolume=2
                        deathcount+=1  
                time = (mcs-relaxtime)/float(10) 
                P_apo.append([cell.id,cell.type,time,dens,deathcount])  
                
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\P_apo-%s'%datetime.now().strftime('%H-%M-%m-%d'),P_apo)
        return
    
    
class neighbourdata(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)   
        
    def start(self):
        self.pW=self.addNewPlotWindow(
        _title='CellCount v Time',
        _xAxisTitle='Time(Frames)',
        _yAxisTitle='Cell Count',
        _xScaleType='linear',
        _yScaleType='log')        
        self.pW.addPlot('cellcount WT',_style='Dots',_color='green',_size=5)   
    def step(self,mcs):
        cellcount_wt=0
        for cell in self.cellList:
            cellcount_wt+=1    
        time = (mcs-relaxtime)/float(10)       
        self.pW.addDataPoint('cellcount WT',time,cellcount_wt)  
        self.pW.savePlotAsData('Cell_Count.txt')
#         self.pW.addDataPoint("cellvol  WT",time,cellvolume/cellcount) 
    def finish(self):
        return     
        
        
        #Tracking
class tracking(SteppableBasePy):   
    
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)     
    def start(self): 
        global trackingfile
        trackingfile=[]
    def step(self,mcs): 
        cellcount=0
        cellvolume=0
        for cell in self.cellList:
            if cell.type==5 or cell.type==6 :
                state=1
            else:
                state=0
            time = (mcs-relaxtime)/float(10)      
            ar1= [cell.xCOM,cell.yCOM,int(time),int(cell.id),int(cell.type),state]
            trackingfile.append(ar1) 
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\\tracking-%s' %datetime.now().strftime('%H-%M-%m-%d'), trackingfile)
        return
