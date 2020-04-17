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
            
            
global minvol
global maxvol
global relaxtime
global adderlist
growthratewt=6.0
growthratescrb=3.2
stiffness_kd=0.9
stiffness_wt=2.0
CI=0.01
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
        volList = []
        self.pW=self.addNewPlotWindow(_title='initial seeding volumes' ,_xAxisTitle='initial volume',_yAxisTitle='N', _xScaleType='linear',_yScaleType='linear')
        self.pW.addHistogramPlot(_plotName='initialVol',_color='green',_alpha=100)# _alpha is transparency 0 is transparent, 255 is opaque        
        for cell in self.cellList:
            if cell.type==2 and 0.95>random.random():
                cell.targetVolume=0
                cell.lambdaVolume=5
#                 cell.targetVolume=random.randint(minvol-100,maxvol+100)    
#                 cell.targetVolume=random.normalvariate(2100, 300)
            else:
                cell.targetVolume=random.normalvariate(1900, 500)
                round(cell.targetVolume)
            if cell.type==1:
                cell.lambdaVolume=stiffness_wt
            elif cell.type==2:
                cell.lambdaVolume=stiffness_kd
            id=cell.id   
            add1=[id, 1900,cell.type]
            adderlist.append(add1)
            volList.append(cell.targetVolume)  
#         self.pW.addHistogram(plot_name='initialVol', value_array=volList, number_of_bins=100)     
#         self.pW.savePlotAsData('initialVol.txt')    












#Growth of cells for mitosis, with noise, also change lambda with time so size reduced.
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
                if cell.type==1 or cell.type==3 or cell.type==5 or cell.type==7:
                    cell.targetVolume+=round(random.normalvariate(growthratewt, 2.5))*np.exp(-(CI/mcs)*((cell.volume-cell.targetVolume)**2))  # a = softness parameter
                    tvadded+=round(random.normalvariate(growthratewt, 2.5))*np.exp(-(CI/mcs)*((cell.volume-cell.targetVolume)**2))
                elif cell.type==2 or cell.type==4 or cell.type==6 or cell.type==8:
                    cell.targetVolume+=round(random.normalvariate(growthratescrb, 2.5))*np.exp(-(CI/mcs)*((cell.volume-cell.targetVolume)**2))  # a = softness parameter
                    tvadded+=round(random.normalvariate(growthratescrb, 2.5))*np.exp(-(CI/mcs)*((cell.volume-cell.targetVolume)**2))
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
                    if x[0]==cell.id and cell.volume-x[1] > (1800-(mcs/60)) and (cell.type!=9 and cell.type!=10):
                        cells_to_divide.append(cell)  
                        Volumeanalysis.append([x[0],(mcs-relaxtime)/float(10),x[1],cell.volume,cell.type])                       
        for cell in cells_to_divide:
            # to change mitosis mode leave one of the below lines uncommented
#             self.divideCellRandomOrientation(cell)
            # self.divideCellOrientationVectorBased(cell,1,0,0)                 # this is a valid option
            # self.divideCellAlongMajorAxis(cell)                               # this is a valid option
            self.divideCellAlongMinorAxis(cell)                               # this is a valid option
 
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
            self.childCell.type=3
            self.parentCell.type=3
        elif self.parentCell.type==2:
            self.childCell.type=4
            self.parentCell.type=4
        elif self.parentCell.type==3:
            self.childCell.type=5   
            self.parentCell.type=5
        elif self.parentCell.type==4:
            self.childCell.type=6
            self.parentCell.type=6
        elif self.parentCell.type==5:
            self.childCell.type=7   
            self.parentCell.type=7
        elif self.parentCell.type==6:
            self.childCell.type=8   
            self.parentCell.type=8  
        elif self.parentCell.type==7:
            self.childCell.type=1  
            self.parentCell.type=1  
        elif self.parentCell.type==8:
            self.childCell.type=2   
            self.parentCell.type=2   
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\parentlist-%s'%datetime.now().strftime('%H-%M-%m-%d'), parentlist)   
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\Volumeanalysis-%s'%datetime.now().strftime('%H-%M-%m-%d'),Volumeanalysis)  
        
        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.cloneAttributes(sourceCell=self.parentCell, targetCell=self.childCell, no_clone_key_dict_list = [attrib1, attrib2] )
        
        
     
        














#Motility
class CellMotilitySteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):        
        # iterating over all cells in simulation        
        for cell in self.cellList:
            break 
            # Make sure ExternalPotential plugin is loaded
            # negative lambdaVecX makes force point in the positive direction
            cell.lambdaVecX=10.1*random.uniform(-1.0,1.0) # force component pointing along X axis 
            cell.lambdaVecY=10.1*random.uniform(-1.0,1.0) # force component pointing along Y axis 
#         cell.lambdaVecZ=0.0 # force component pointing along Z axis 
        
        
    def step(self,mcs):
        
        for cell in self.cellList:
            
            cell.lambdaVecX=10.1*random.uniform(-1.0,1.0) # force component pointing along X axis 
            cell.lambdaVecY=10.1*random.uniform(-1.0,1.0) # force component pointing along Y axis 
            # print 'cell.lambdaVecX=',cell.lambdaVecX,' cell.lambdaVecY=',cell.lambdaVecY




class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
       global LocalDensity
       global P_apo
       global cumulitiveapop
       cumulitiveapop=[]
       LocalDensity=[]
       P_apo=[]
#        self.pW=self.addNewPlotWindow(
#        _title=' cumulitive apoptosis',
#        _xAxisTitle='Time(Frames)',
#        _yAxisTitle='apopstosis',
#        _xScaleType='linear',
#        _yScaleType='linear')
        
#        self.pW.addPlot('apoptosis WT',_style='Dots',_color='green',_size=5)   
#        self.pW.addPlot('aposptosis scrb',_style='Dots',_color='red',_size=5) 
       
    def step(self,mcs):
        apopwt=0
        apopscrb=0
        deathcount=0
        dens=0
        totdens=0
        cellcount=0     
        if mcs>relaxtime:
            for cell in self.cellList: 
                if cell.type!=9 and cell.type!=10:
#                     if (cell.xCOM>200 or cell.xCOM<1000) and (cell.yCOM>200 or cell.xCOM<1400):
                    totdens=0
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
                    if cell.type==2 or cell.type==4 or cell.type==6 or cell.type==8:
                        if 0.00072194/(1 + np.exp(-509.4*(dens*3- 0.0067)))>random.random():#density dependent death
                        
                            cell.type=10
                            cell.targetVolume=0
                            cell.lambdaVolume=2
                            deathcount+=1
                            apopscrb+=1
                    if cell.type==1 or cell.type==3 or cell.type==5 or cell.type==7:
                        if 0.00033074/(1 + np.exp(-235.8*(dens*3- (0.0152))))>random.random():#density dependent death
                        
                            cell.type=9
                            cell.targetVolume=0
                            cell.lambdaVolume=2
                            deathcount+=1
                            apopwt+=1
                time = (mcs-relaxtime)/float(10) 
                if cell.type==2 or cell.type==4 or cell.type==6 or cell.type==8:           
                    LocalDensity.append([cell.id,cell.type,time,totdens])
                elif cell.type==1 or cell.type==3 or cell.type==5 or cell.type==7:  
                    LocalDensity.append([cell.id,cell.type,time,totdens])  
                P_apo.append([cell.id,cell.type,time,totdens,deathcount])  
                cumulitiveapop.append([time,apopscrb,apopwt])  

                
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\LocalDensity-%s'%datetime.now().strftime('%H-%M-%m-%d'), LocalDensity)   
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\P_apo-%s'%datetime.now().strftime('%H-%M-%m-%d'),P_apo)
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\cumulitiveapop-%s'%datetime.now().strftime('%H-%M-%m-%d'),cumulitiveapop)
        return



        
#  Apoptosis
class DeathSteppablePerimiter(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
       global perimeterarray
       global P_apo
       global cumulitiveapop
       cumulitiveapop=[]
       perimeterarray=[]
       P_apo=[]
#        self.pW=self.addNewPlotWindow(
#        _title=' cumulitive apoptosis',
#        _xAxisTitle='Time(Frames)',
#        _yAxisTitle='apopstosis',
#        _xScaleType='linear',
#        _yScaleType='linear')
        
#        self.pW.addPlot('apoptosis WT',_style='Dots',_color='green',_size=5)   
#        self.pW.addPlot('aposptosis scrb',_style='Dots',_color='red',_size=5) 
       
    def step(self,mcs):
        apopwt=0
        apopscrb=0
        deathcount=0
        cellcount=0 
        typecell=0
        if mcs>relaxtime:
            for cell in self.cellList: 
                if cell.type==2 or cell.type==4 or cell.type==6 or cell.type==8:
                    typecell=2
                else:
                    typecell=1
                neighbourpixel=[]
                perimeterpercentage=0
                if cell.type!=9 and cell.type!=10:
#                     if (cell.xCOM>200 or cell.xCOM<1000) and (cell.yCOM>200 or cell.xCOM<1400):
#                     print "cell ID = ", cell.id, "cell Type = ", cell.type
                    for pixel in self.getCopyOfCellBoundaryPixels(cell):
                        a=self.point3DToNumpy(pixel)
                        centrecell = self.cellField[a[0],a[1],a[2]] #x, y, z are the pixel coordinates 
#                         print "cell ID = ", centrecell.id, "cell Type = ", centrecell.type
# #                        
#                         for neigbouringpixel in self.getPixelNeighborsBasedOnNeighborOrder(pixel,1):
#                             print neigbouringpixel
# #                             aa=self.point3DToNumpy(neigbouringpixel)
# #                             print aa
                        
                        
                        
                        leftpixel = self.cellField[a[0]-1,a[1],a[2]]
                        if hasattr(leftpixel, "id"):
                            if centrecell.id!=leftpixel.id:
                                neighbourpixel.append(leftpixel.type)
                        else: 
                            neighbourpixel.append(0)     
                        rightpixel = self.cellField[a[0]+1,a[1],a[2]]
                        if hasattr(rightpixel, "id"):
                            if centrecell.id!=rightpixel.id:
                                neighbourpixel.append(rightpixel.type)
                        else:   
                            neighbourpixel.append(0)  
                        uppixel = self.cellField[a[0],a[1]+1,a[2]]
                        if hasattr(uppixel, "id"):
                            if centrecell.id!=uppixel.id:
                                neighbourpixel.append(uppixel.type)
                        else:   
                            neighbourpixel.append(0) 
                        downpixel = self.cellField[a[0],a[1]-1,a[2]]
                        if hasattr(downpixel, "id"):
                            if centrecell.id!=downpixel.id:
                                neighbourpixel.append(downpixel.type) 
                        else:   
                            neighbourpixel.append(0)                                
            
#                 print "cell ",cell.id," length of perimeter = ",len(neighbourpixel),"length of shared perimiter being same type =", neighbourpixel.count(cell.type)                 
                
                if typecell==2 and len(neighbourpixel)>0 :
                    
                    perimeterpercentage=((neighbourpixel.count(1)+ neighbourpixel.count(3) + neighbourpixel.count(5) + neighbourpixel.count(7)) / len(neighbourpixel))                
                    if 3*(0.000416*perimeterpercentage) > random.random():
                       
                        cell.type=10
                        cell.targetVolume=0
                        cell.lambdaVolume=2
                        deathcount+=1
                        apopscrb+=1
                
                elif typecell==1 and len(neighbourpixel)>0 :
                    perimeterpercentage=((neighbourpixel.count(2)+neighbourpixel.count(4)+neighbourpixel.count(6)+neighbourpixel.count(8))/len(neighbourpixel))
#winner-loser shared contact dependent death                    
                    if 10*(0.0000416*perimeterpercentage)>random.random(): 
                        
                        
                            cell.type=9
                            cell.targetVolume=0
                            cell.lambdaVolume=2
                            deathcount+=1
                            apopwt+=1
                time = (mcs-relaxtime)/float(10) 
                P_apo.append([cell.id,cell.type,time,perimeterpercentage,deathcount])  
            time = (mcs-relaxtime)/float(10) 
            perimeterarray.append([cell.id,cell.type,time,perimeterpercentage])   
            cumulitiveapop.append([time,apopscrb,apopwt])  

                
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\perimeterarray-%s'%datetime.now().strftime('%H-%M-%m-%d'), perimeterarray)   
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\P_apo-%s'%datetime.now().strftime('%H-%M-%m-%d'),P_apo)
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\cumulitiveapop-%s'%datetime.now().strftime('%H-%M-%m-%d'),cumulitiveapop)
        return       
        
        








# #Apoptosis
# class DeathSteppable(SteppableBasePy):
#     def __init__(self,_simulator,_frequency=5):
#         SteppableBasePy.__init__(self,_simulator,_frequency)
#     def start(self):
#        global LocalDensity
#        global P_apo
#        LocalDensity=[]
#        P_apo=[]
#     def step(self,mcs):
#         deathcount=0
#         dens=0
#         totdens=0
#         cellcount=0     
#         if mcs>relaxtime:
#             for cell in self.cellList: 
#                 if cell.type!=9 and cell.type!=10:
# #                     if (cell.xCOM>200 or cell.xCOM<1000) and (cell.yCOM>200 or cell.xCOM<1400):
#                     totdens=0
#                     cellneighbours=[]
#                     nieghbourvolume=[]
#                     numberOfCellNeighbours=[]
#                     dens=(1/cell.volume)
#                     for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
#                         if neighbor:
#                             cellneighbours.append(neighbor.id)   
#                             nieghbourvolume.append(neighbor.volume)
#                             dens+=(1/neighbor.volume)
#                     cn= (str(w) for w in cellneighbours)
#                     den= (str(k) for k in nieghbourvolume)
#                     li=list(nieghbourvolume)
#                     ne=len(cellneighbours)
#                     totdens+=dens
#                     cellcount+=1
#                     if cell.type==2 or cell.type==4 or cell.type==6 or cell.type==8:
#                         if 0.00072194/(1 + np.exp(-509.4*(dens*9- 0.0067)))>random.random():#density dependent death
                        
#                             cell.type=10
#                             cell.targetVolume=0
#                             cell.lambdaVolume=2
#                             deathcount+=1
#                     elif cell.type==1 or cell.type==3 or cell.type==5 or cell.type==7:
#                         if 0.00033074/(1 + np.exp(-235.8*(dens*9- 0.0152)))>random.random():#density dependent death
                        
#                             cell.type=9
#                             cell.targetVolume=0
#                             cell.lambdaVolume=2
#                             deathcount+=1
#                 time = (mcs-relaxtime)/float(10) 
#                 if cell.type==2 or cell.type==4 or cell.type==6 or cell.type==8:           
#                     LocalDensity.append([cell.id,cell.type,time,totdens])
#                 elif cell.type==1 or cell.type==3 or cell.type==5 or cell.type==7:  
#                     LocalDensity.append([cell.id,cell.type,time,totdens])  
#                 P_apo.append([cell.id,cell.type,time,totdens,deathcount])  
                
#     def finish(self):
#         np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\LocalDensity-%s'%datetime.now().strftime('%H-%M-%m-%d'), LocalDensity)   
#         np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\P_apo-%s'%datetime.now().strftime('%H-%M-%m-%d'),P_apo)
#         return
        
        
        
        
 

# class DeathSteppable(SteppableBasePy):
#     def __init__(self,_simulator,_frequency=5):
#         SteppableBasePy.__init__(self,_simulator,_frequency)
#     def start(self):
#        global P_apo
#        global LocalDensity
#        LocalDensity=[]
#        P_apo=[]
#     def step(self,mcs):
#         deathcount=0
#         dens=0
#         totdens=0
#         cellcount=0     
#         for cell in self.cellList: 
#             if cell.type!=9 and cell.type!=10:
#                 dens=(1/cell.volume)
#                 if cell.type==2 or cell.type==4 or cell.type==6 or cell.type==8:
#                     if 0.00072194/(1 + np.exp(-509.4*(dens*9- 0.0026)))>random.random():#density dependent death
                    
#                         cell.type=10
#                         cell.targetVolume=0
#                         cell.lambdaVolume=2
#                         deathcount+=1
#                 elif cell.type==1 or cell.type==3 or cell.type==5 or cell.type==7:
#                     if 0.00033074/(1 + np.exp(-235.8*(dens*9- 0.0152)))>random.random():#density dependent death
                    
#                         cell.type=9
#                         cell.targetVolume=0
#                         cell.lambdaVolume=2
#                         deathcount+=1
#             time = (mcs-relaxtime)/float(10) 
#             if cell.type==2 or cell.type==4 or cell.type==6 or cell.type==8:           
#                 LocalDensity.append([cell.id,cell.type,time,dens])
#             elif cell.type==1 or cell.type==3 or cell.type==5 or cell.type==7:  
#                 LocalDensity.append([cell.id,cell.type,time,dens])  
#             P_apo.append([cell.id,cell.type,time,totdens,deathcount])  
            
#     def finish(self):
#         np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\LocalDensity-%s'%datetime.now().strftime('%H-%M-%m-%d'), LocalDensity)   
#         np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\P_apo-%s'%datetime.now().strftime('%H-%M-%m-%d'),P_apo)
#         return
# #Apoptosis
# class DeathSteppable_Elasticity(SteppableBasePy):
#     def __init__(self,_simulator,_frequency=10):
#         SteppableBasePy.__init__(self,_simulator,_frequency)
#     def start(self):
#        global P_apo
#        P_apo=[]
#     def step(self,mcs):
#         deathcount=0
#         dens=0
#         totdens=0
#         cellcount=0
        
#         if mcs>relaxtime:
#             for cell in self.cellList: 
#                 if cell.type!=9 and cell.volume>1200:
#                         cellcount+=1
#                         if cell.type==2 or cell.type==4 or cell.type==6 or cell.type==8:
#                             if np.exp(-0.5*((cell.volume-cell.targetVolume))) < random.random()/10000:#density dependent death
#                                 cell.type=9
#                                 cell.targetVolume=0
#                                 cell.lambdaVolume=2
#                                 deathcount+=1
#                 time = (mcs-relaxtime)/float(10)                
#                 P_apo.append([cell.id,cell.type,time,deathcount])  
#     def finish(self):   
#         np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\P_apo',P_apo)
#         return
        
            
















#CellCount
class neighbourdata(SteppableBasePy):
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
        self.pW.addPlot('cellcount scrb',_style='Dots',_color='red',_size=5) 
#         self.pW=self.addNewPlotWindow(_title='Cellvolume v Time',_xAxisTitle='Time(Frames)',_yAxisTitle='Cell volume',_xScaleType='linear',_yScaleType='log')
#         self.pW.addPlot('cellvol WT',_style='Dots',_color='blue',_size=5) 
    def step(self,mcs):
        cellcount_wt=0
        cellcount_scrb=0
        for cell in self.cellList:
            if cell.type==1 or cell.type==3 or cell.type==5 or cell.type==7:
                cellcount_wt+=1
            elif cell.type==2 or cell.type==4 or cell.type==6 or cell.type==8: 
                cellcount_scrb+=1    
        time = (mcs-relaxtime)/float(10)       
        self.pW.addDataPoint('cellcount WT',time,cellcount_wt)  
        self.pW.savePlotAsData('Cell_Count.txt')
        self.pW.addDataPoint('cellcount scrb',time,cellcount_scrb)  
        self.pW.savePlotAsData('Cell_Count_scrb.txt')
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
            if cell.type==9 or cell.type==10:
                state=1
            else:
                state=0
            time = (mcs-relaxtime)/float(10)      
            ar1= [cell.xCOM,cell.yCOM,int(time),int(cell.id),int(cell.type),state]
            trackingfile.append(ar1) 
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\\tracking-%s' %datetime.now().strftime('%H-%M-%m-%d'), trackingfile)
        return


 








# class P_locdens(SteppableBasePy):
#     def __init__(self,_simulator,_frequency=10):
#         SteppableBasePy.__init__(self,_simulator,_frequency)
        
#     def start(self):
#         self.pW=self.addNewPlotWindow(_title='P_div',_xAxisTitle='Avg Local Density',_yAxisTitle='P_div',_xScaleType='linear',_yScaleType='linear')
#         self.pW.addPlot('p_div_wt',_style='Dots',_color='blue',_size=5)   

#     def step(self,mcs):
#         p_div=0
#         totdens=0
#         cellcount=0
#         dens=0
#         Div_list=len(parentlist)
#         if mcs>relaxtime:    
#             for cell in self.cellList:   
#                 if (cell.xCOM>200 or cell.xCOM<1000) and (cell.yCOM>200 or cell.xCOM<1400):
#                     cellneighbours=[]
#                     nieghbourvolume=[]
#                     numberOfCellNeighbours=[]
#                     dens=(1/cell.volume)
#                     for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
#                         if neighbor:
#                             cellneighbours.append(neighbor.id)   
#                             nieghbourvolume.append(neighbor.volume)
#                             dens+=(1/neighbor.volume)
#                     li=list(nieghbourvolume)
#                     totdens+=dens
#                     cellcount+=1             
#             self.pW.addDataPoint('p_div_wt',totdens/cellcount,Div_list/cellcount)  
#             self.pW.savePlotAsData('p_div_wt.txt')   
            
            


class cleanup(SteppableBasePy):   
    
    
    def __init__(self,_simulator,_frequency=100):
        SteppableBasePy.__init__(self,_simulator,_frequency)     
        global extrusions
        extrusions=[]
    def step(self,mcs): 
        if mcs>relaxtime:
            for cell in self.cellList:
                if cell.volume<500:
                    extrusions.append([cell.id,cell.type,(mcs-relaxtime)/float(10)]) 
                if cell.type==9 or cell.type==10 or cell.volume<500 :
                    self.deleteCell(cell)
                    
                
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\extrusions-%s'%datetime.now().strftime('%H-%M-%m-%d'), extrusions)   
        return       
  #Density Dependent death.
