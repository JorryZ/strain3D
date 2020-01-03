# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 13:16:10 2019

@author: bieZY (e0348827@u.nus.edu)
Objective: calculate the 3D strain based on the b-Spline Fourier model
Method: 
Modules: numpy, scipy, motionSegmentation, trimesh

History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: jorry.zhengyu@gmail.com         20NOV2019           -V1.0.0 Created, test version
  Author: jorry.zhengyu@gmail.com         07Dec2019           -V1.0.1 function strain3D.centerPointFit: add dimDist, ajudgment for inner and outer points
  Author: jorry.zhengyu@gmail.com         03JAN2020           -V1.0.2 function strain3D.centerPointFit: add more option for stlSample
"""
print('strainMyocardium test version 1.0.2')

import sys
import numpy as np
import scipy.spatial as spatial
from scipy.spatial.distance import cdist
import medImgProc
import motionSegmentation.BsplineFourier as BsplineFourier
import trimesh

def normalize_v3(arr):
    ''' Normalize a numpy array of 3 component vectors shape=(n,3) '''
    lens = np.sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens                
    return arr

class strain3D:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self):
        '''
        Initialize motionConstrain class
        '''
        self.time=None
        self.origin=None
        self.spacing=None   #grid spacing (x,y,z)
        self.period=None
        self.coefMat=[]
        self.shape=None
        self.coordMat=None
        self.bsFourierPath=None
        
        self.dimlen=None    #physical length of each (x,y,z) pixel
        self.sampleCoord =None
        self.sampleNormal=None
        self.bsFourier=BsplineFourier.BsplineFourier()
        
        self.minDist=None   #the min distances among points
        self.minAxes=None   #the min distances from each axes
        self.assisCoord=None
        self.assisCoordX=None
        self.assisCoordY=None
        self.assisCoordZ=None
        self.vector=None
        
        self.coordCoef=[]
        self.coordMotion=None
        self.tensor=None
        self.strain=None
        
        self.stlVertice = None
        self.faceNormal = None
        self.faceCenter = None
        self.innerFaceNormal = None
        self.innerFaceCenter = None
        self.outerFaceNormal = None
        self.outerFaceCenter = None
        self.centerPoint = None
        self.centerLine = None
        self.longAxis = None
        self.longAxisMin = None
        self.longAxisMax = None
        
        self.longitAxis = None
        self.circumAxis = None
        self.radialAxis = None
        
        self.longitStrain = None
        self.circumStrain = None
        self.radialStrain = None
        self.eigenVector = None
        self.eigenValue = None
        print('strain3D_init_ done')
        
    def initialize(self,coefFile=None,origin=None,shape=None,spacing=None,period=None,time=None,fourierFormat=None,fileScale=1.,delimiter=' ',getSampleCoord=False,spacingDivision=[1.,1.5],gap=0):   
        if type(coefFile)!=type(None):
            self.bsFourierPath = coefFile
            self.bsFourier.readFile(coefFile=coefFile,origin=origin,shape=shape,spacing=spacing,fourierFormat=fourierFormat,delimiter=delimiter)
            self.origin=self.bsFourier.origin
            self.spacing=self.bsFourier.spacing[:3]
            self.period=self.bsFourier.spacing[3]
            self.coefMat=self.bsFourier.coef.copy()
            self.shape=self.coefMat.shape
            self.coordMat = np.mgrid[self.origin[0]:(self.origin[0]+(self.shape[0]-0.1)*self.spacing[0]):self.spacing[0], self.origin[1]:(self.origin[1]+(self.shape[1]-0.1)*self.spacing[1]):self.spacing[1],self.origin[2]:(self.origin[2]+(self.shape[2]-0.1)*self.spacing[2]):self.spacing[2]].reshape(3,*self.shape[:3]).transpose(1,2,3,0)
            #self.coefMat=self.bsFourier.coef.reshape(-1, order='F')
        if getSampleCoord==True:
            self.sampleCoord=np.array(self.bsFourier.samplePoints(spacingDivision=spacingDivision[0],gap=gap))  #self.bsFourier.samplePoints()
            print('Initialization: sampleCoord has %d points.'%(len(self.sampleCoord)))
        
        if type(time)!=type(None):
            self.time=time
            print('Initialization: time=',self.time)
        
        if type(self.bsFourier.coef)==np.ndarray:
          print('Initialization: shape=',self.bsFourier.coef.shape)
          print('Initialization: spacing=',self.bsFourier.spacing)
          print('Initialization: origin=',self.bsFourier.origin)

    def distCalc(self,coordA=None, coordB=None, mode='min'):
        '''
        calculate the distance among coordA: min, max (to do)
        calculate the distance between coordA and coordB: min (to do), max (to do)
        '''
        # need test
        if type(coordB)==type(None) and mode=='min':
            if type(coordA)==type(None):
                coordA=self.sampleCoord
            flag = 0
            #p1=np.zeros((1,3))
            temp=coordA[0]
            p1=np.reshape(temp,(1,np.product(temp.shape)))
            p2=coordA
            distList = spatial.distance.cdist(p1,p2)
            temp = np.reshape(distList,(np.product(distList.shape),))
            minDistValue = np.min(temp[np.nonzero(temp)])
            minDistIndex = np.where(temp==minDistValue)
            temp=minDistIndex[0]
            minDistIndex = temp[0]
            '''
            minDistValue = heapq.nsmallest(2, temp)[-1]
            minDistIndex = np.where(temp==minDistValue)
            minDistIndex = minDistIndex[0]
            '''
            #minDistValue = np.min(temp)     #find the min value
            #minDistIndex = np.argmin(temp)  #find the index of the min value
            
            while flag!=1:
                temp=coordA[minDistIndex]
                p1=np.reshape(temp,(1,np.product(p1.shape)))
                #p2=np.delete(self.sampleCoord,minDistIndex,axis=0)
                distList = spatial.distance.cdist(p1,p2)
                temp = np.reshape(distList,(np.product(distList.shape),))
                minTemp = np.min(temp[np.nonzero(temp)])
                if minTemp < minDistValue:
                    minDistValue = minTemp
                    minDistIndex = np.where(temp==minDistValue)
                    temp=minDistIndex[0]
                    minDistIndex = temp[0]
                else:
                    flag=1
            return minDistValue
        
        if type(coordB)!=type(None):
            temp = spatial.distance_matrix(coordA,coordB)
            temp = temp.transpose()
            dist=temp.tolist()
            distValue=[]
            distIndex=[]
            if mode=='min':
                for i in range(len(coordB)):
                    distValue.append(min(dist[i]))
                    distIndex.append(dist[i].index(distValue[i]))
            elif mode=='max':
                for i in range(len(coordB)):
                    distValue.append(max(dist[i]))
                    distIndex.append(dist[i].index(distValue[i]))
            return distValue,distIndex
        
    def minAxesCalc(self):
        '''
        calculate the min distances from each axes (x,y,z), 
        '''
        self.minAxes=np.zeros(self.spacing.shape)
        for axes in range(len(self.spacing)):
            temp=self.sampleCoord[:,axes]
            data=np.reshape(temp,(np.product(temp.shape),1))
            self.minAxes[axes]=self.distCalc(coordA=data,mode='min')
            
    def assisCoordCalc(self, division=5.):
        '''
        get assisCoord based on the minAxes, for calculation of displacement
        '''
        assisCoord = self.sampleCoord
        dim = len(self.spacing)
        for i in range(dim-1):
            assisCoord=np.concatenate((assisCoord,self.sampleCoord),axis=1)
        if type(self.minAxes)!=type(None):
            minValue = self.minAxes
            for axes in range(dim):
                assisCoord[:,axes*dim+axes] = assisCoord[:,axes*dim+axes]-minValue[axes]/division
            self.assisCoord=assisCoord.copy()
        else:
            minValue = self.minDist
            for axes in range(dim):
                assisCoord[:,axes*dim+axes] = assisCoord[:,axes*dim+axes]-minValue/division
            self.assisCoord=assisCoord.copy()
            self.assisCoordX = self.assisCoord[:,0:dim]
            self.assisCoordY = self.assisCoord[:,dim:2*dim]
            self.assisCoordZ = self.assisCoord[:,2*dim:3*dim]
            
        
    def coordMotionCalc(self,sampleCoord=None,time=None):
        '''
        calculate the point movement over the cardiac cycle
        parameter: time is the referance time in the cardiac cycle
        '''
        if type(sampleCoord)==type(None):
            dim = len(self.spacing)
            temp=np.concatenate((self.sampleCoord,self.assisCoord),axis=1)
            sampleCoord=np.reshape(temp,(len(temp)*(dim+1),dim),order='C')
            
        if type(time)==type(None):
            time=np.arange(self.period)
        
        coordCoef=[]
        for m in range(len(sampleCoord)):
            coef=self.bsFourier.getRefCoef(sampleCoord[m])
            coordCoef.append(coef.copy())
        self.coordCoef=np.array(coordCoef.copy())   #[sampleCoord, ft, uvw]
        P=self.period
        for n in time:  #time points to calculate coordinates
            temp=sampleCoord.copy()
            temp+=self.coordCoef[:,0]
            for m in range(int(len(self.coordCoef[0])/2)):
                #pointsCoef[:,0] is 0
                
                temp+=self.coordCoef[:,m+1]*np.cos((m+1.)*2.*np.pi/P*n)+self.coordCoef[:,int(len(self.coordCoef[0])/2)+m+1]*np.sin((m+1.)*2.*np.pi/P*n)
            if type(self.coordMotion)==type(None):
                self.coordMotion=temp
            else:
                self.coordMotion=np.concatenate((self.coordMotion,temp),axis=1)
        
    def vectorCalc(self):
        '''
        calculate displacement AX, AY, AZ
        '''
        vector=[]
        dim = len(self.spacing)
        temp=self.coordMotion
        size=self.coordMotion.shape
        #coordMotion=np.reshape(temp,(int(size[0]/(dim+1)),int(size[1]*(dim+1))),order='C')
        length = int(size[0]/(dim+1))
        coordMotionA = temp[0:1*length,:]
        coordMotionX = temp[1*length:2*length,:]
        coordMotionY = temp[2*length:3*length,:]
        coordMotionZ = temp[3*length:4*length,:]
        
        vectorX=coordMotionA-coordMotionX
        vectorY=coordMotionA-coordMotionY
        vectorZ=coordMotionA-coordMotionZ
        vector.append(vectorX)
        vector.append(vectorY)
        vector.append(vectorZ)
        self.vector=vector.copy()
        #self.vector=np.array(vector.copy())
        
    def stlNormal(self, stlName):
        # get the normal of the stl surface
        stlData=trimesh.load(stlName)
        vertices = stlData.vertices
        faces = stlData.faces
        norm = np.zeros( vertices.shape, dtype=vertices.dtype )
        tris = vertices[faces]
        n = np.cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )    #normal of each triangle
        normalize_v3(n)
        #n=n*(-1)
        faceCenter = np.mean(tris,axis=1)
        self.stlVertice = vertices.copy()
        self.faceNormal = n.copy()
        self.faceCenter = faceCenter.copy()
        
        #normal of each vertice
        norm[ faces[:,0] ] += n
        norm[ faces[:,1] ] += n
        norm[ faces[:,2] ] += n
        normalize_v3(norm)
    
    def centerPointFit(self, dataName=None, dimlen=None, longAxis='z', badSlice=[5,3], stlSample=False):
        '''
        fit the center point of circle-like points
        longAxis is the direction of center line
        badSlice: abandoned slice to fit center, apex part first, basal part second
            seperate inner surface points and outer surface points
        stlSample: True: assign sampleData to sampleCoord (vertices), 'cut': sampleData and abandon part of apex and basal (vertices); False: no;
                    'innercut': innerFaceCenter and abandon part of apex and basal (faces); 'outercut': outerFaceCenter and abandon part of apex and basal (faces)
        '''
        if type(dataName)==type(None):
            try:
                sampleData = self.stlVertice.copy()
            except:
                print('Warning: Function centerFit - dataName is None and sampleData is None')
                sys.exit()
        
        elif dataName[-3:]=='vtk':
            if type(dimlen)==type(None):
                print('Warning: Function centerFit - source of data from VTK file, please specify a VTK file and a dimlen (dimension length)!')
                sys.exit()
            sampleData=medImgProc.imread(dataName)
            sampleData.dim=['z','y','x']
            sampleData.dimlen=dimlen
            sampleData.rearrangeDim(['x','y','z'])
            
        elif dataName[-3:]=='stl':
            if type(dimlen)==type(None):
                print('Warning: Function centerPointFit - source of data from STL file, please specify a STL file and a dimlen (dimension length)!')
                sys.exit()
            stlData=trimesh.load(dataName)
            sampleData=stlData.vertices
        if stlSample==True:
            self.sampleCoord=sampleData.copy()
        
        dimDist = {'x':0,'y':1,'z':2,'X':0,'Y':1,'Z':2}
        dim = dimDist[longAxis]
        
        faceCenter=self.faceCenter.copy()
        faceNormal=self.faceNormal.copy()
        longAxisMin = np.min(sampleData[:,dim])
        longAxisMax = np.max(sampleData[:,dim])
        sliceLoc = longAxisMin + badSlice[0]*dimlen[dim]            #location of start slice
        centerPoint = []
        innerFaceCenter = None
        innerFaceNormal = None
        outerFaceCenter = None
        outerFaceNormal = None
        while sliceLoc < (longAxisMax-badSlice[1]*dimlen[dim]):       #location of end slice
            #calculate center point, and seperate inner and outer face center
            points = np.array([sampleData[i,:] for i in range(len(sampleData)) if sampleData[i,dim]<(sliceLoc+dimlen[dim]) and sampleData[i,dim]>(sliceLoc-dimlen[dim])])
            center = np.mean(points,axis=0)
            centerPoint.append(center)
            sliceLoc=sliceLoc+2*dimlen[dim]
            
            tempCenter = center.reshape((1,3))
            index = np.array([i for i in range(len(faceCenter)) if faceCenter[i,dim]<(sliceLoc+dimlen[dim]) and faceCenter[i,dim]>(sliceLoc-dimlen[dim])])
            temp1=faceCenter[index]
            dist1=spatial.distance_matrix(temp1,tempCenter)
            temp2=temp1-faceNormal[index]
            dist2=spatial.distance_matrix(temp2,tempCenter)
            difDist=np.array(dist1-dist2)
            inner = np.array([i for i in range(len(difDist)) if difDist[i]<0])
            outer = np.array([i for i in range(len(difDist)) if difDist[i]>0])
            if len(inner)!=0:
                try:
                    innerFaceCenter=np.concatenate((innerFaceCenter,faceCenter[index[inner]]),axis=0)
                    innerFaceNormal=np.concatenate((innerFaceNormal,faceNormal[index[inner]]),axis=0)
                except:
                    innerFaceCenter=np.array(faceCenter[index[inner]])
                    innerFaceNormal=np.array(faceNormal[index[inner]])
            if len(outer)!=0:
                try:
                    outerFaceCenter=np.concatenate((outerFaceCenter,faceCenter[index[outer]]),axis=0)
                    outerFaceNormal=np.concatenate((outerFaceNormal,faceNormal[index[outer]]),axis=0)
                except:
                    outerFaceCenter=np.array(faceCenter[index[outer]])
                    outerFaceNormal=np.array(faceNormal[index[outer]])
            if stlSample=='cut':
                if type(self.sampleCoord)==type(None):
                    self.sampleCoord=points
                else:
                    self.sampleCoord=np.concatenate((self.sampleCoord,points),axis=0)
            
        self.centerPoint = np.array(centerPoint.copy())
        self.innerFaceCenter = np.array(innerFaceCenter.copy())
        self.innerFaceNormal = np.array(innerFaceNormal.copy())
        self.outerFaceCenter = np.array(outerFaceCenter.copy())
        self.outerFaceNormal = np.array(outerFaceNormal.copy())
        self.longAxis = dim
        self.longAxisMin = longAxisMin
        self.longAxisMax = longAxisMax
        
        if stlSample=='innercut':
            self.sampleCoord=self.innerFaceCenter
        elif stlSample=='outercut':
            self.sampleCoord=self.outerFaceCenter
        elif stlSample==True:
            self.sampleCoord=sampleData
        
    def centerLineFit(self, data=None):
        '''
        fit center line
        '''
        if type(data)==type(None):
            try:
                data = self.centerPoint.copy()
            except:
                print('Warning: Function centerLineFit - data is None')
                sys.exit()
        datamean = data.mean(axis=0)
        dist=spatial.distance_matrix(data,datamean.reshape((1,3)))
        distMax = np.max(dist)
        # Do an SVD on the mean-centered data.
        u,s,vh = np.linalg.svd(data - datamean)
        
        linepts = vh[0] * np.mgrid[-1*distMax:distMax:2j][:, np.newaxis]
        linepts += datamean
        self.centerLine=linepts.copy()
        
        '''
        import matplotlib.pyplot as plt
        import mpl_toolkits.mplot3d as m3d
        
        ax = m3d.Axes3D(plt.figure())
        ax.scatter3D(*data.T)
        ax.plot3D(*linepts.T)
        plt.show()
        '''
        
    def tensorNumericalCalc(self,sampleCoord=None,time=None,ifMinAxesCalc=False,division=5.):
        '''
        numerical method to calculate tensor
        time needs include at least two time points: reference time, deformed time
        '''
        if type(sampleCoord)==type(None):
            if type(self.sampleCoord)!=type(None):
                sampleCoord=self.sampleCoord
            else:
                print('error: please input sampleCoord to calculate strain')
                sys.exit()
        else:
            self.sampleCoord=sampleCoord
        if type(time)==type(None):
            time=np.arange(self.period)
            
        dim = len(self.spacing)
        tensor = []
        
        self.minDist=self.distCalc(coordA=sampleCoord,mode='min')
        if ifMinAxesCalc:
            self.minAxesCalc()
        self.assisCoordCalc(division=division)
        sampleCoord = np.concatenate((sampleCoord,self.assisCoordX),axis=0)
        sampleCoord = np.concatenate((sampleCoord,self.assisCoordY),axis=0)
        sampleCoord = np.concatenate((sampleCoord,self.assisCoordZ),axis=0)
        self.coordMotionCalc(sampleCoord=sampleCoord,time=time)
        self.vectorCalc()
        vectorX=self.vector[0]
        vectorY=self.vector[1]
        vectorZ=self.vector[2]
        #refTime = time[0]
        #arbTime = time[1]
        disRef = np.zeros((len(vectorX),dim,dim))
        disArb = np.zeros((len(vectorX),dim,dim))
        #disRef = np.array(vectorX[:,0*dim:1*dim],vectorY[:,0*dim:1*dim],vectorZ[:,0*dim:1*dim])   #[sampleCoord, xyz]
        #disArb = np.array(vectorX[:,1*dim:2*dim],vectorY[:,1*dim:2*dim],vectorZ[:,1*dim:2*dim])
        for m in range(len(vectorX)):
            disRef[m,:,0]=vectorX[m,0*dim:1*dim]
            disRef[m,:,1]=vectorY[m,0*dim:1*dim]
            disRef[m,:,2]=vectorZ[m,0*dim:1*dim]
            disArb[m,:,0]=vectorX[m,1*dim:2*dim]
            disArb[m,:,1]=vectorY[m,1*dim:2*dim]
            disArb[m,:,2]=vectorZ[m,1*dim:2*dim]
            
            A = disRef[m]   #.transpose()
            B = disArb[m]   #.transpose()
            # tensor = B/A, matrix right division
            temp = np.linalg.lstsq(A.T, B.T)[0].T
            tensor.append(temp)
        self.tensor = np.array(tensor.copy())
        
    def tensorTheoreticalCalc(self, sampleCoord=None, time=None):
        '''
        theoretical method to calculate tensor, based on BSF model
        time is one deformed time, reference time is fixed by the BSF.txt
        '''
        if type(sampleCoord)==type(None):
            if type(self.sampleCoord)!=type(None):
                sampleCoord=self.sampleCoord
            else:
                print('error: please input sampleCoord to calculate strain')
                sys.exit()
        if type(time)==type(None):
            time=np.arange(self.period)
        #dCList,CIndList=self.bsFourier.getdC(sampleCoord)
        resultBC=[]
        for coord in sampleCoord:
            result=self.bsFourier.getRefCoef(coord)     #fourier coefficient
            resultBC.append(result)
        resultdX=self.bsFourier.getdX(sampleCoord)   # get dUdX, [sampleCoord, dxdydz, FT, dudvdw] matrix, like [sampleCoord, 3, FT, 3]
        tensor = []
        dim = len(self.spacing)
        
        ft = self.shape[-2]
        A_matrix = np.zeros((1,ft))
        for m in range(int(ft/2)):
            A_matrix[0,m+1] = np.cos(time[0]*(m+1)*2*np.pi/self.period)
            A_matrix[0,m+1+int(ft/2)] = np.sin(time[0]*(m+1)*2*np.pi/self.period)
        
        for i in range(len(resultdX)):
            dX=resultdX[i]
            #tensor.append([])
            temp = []
            for m in range(dim):
                #temp.append([])
                B_matrix=dX[m]
                #test
                #B_matrix=dX
                temp.append(A_matrix.dot(B_matrix))
            temp2=np.reshape(np.array(temp),(dim,dim)).transpose()+np.diag((np.ones(dim)))
            tensor.append(temp2.copy())
        self.tensor = np.array(tensor.copy())
        #print(np.max(np.abs(self.tensor[:,1,0])))
        
    def strainCalc(self, tensor=None):
        if type(tensor)==type(None):
            tensor=self.tensor
        strain=[]
        dim = len(self.spacing)
        
        for i in range(len(tensor)):
            temp=tensor[i]
            tempT=temp.transpose()
            I_matrix=np.identity(dim)
            strain.append(0.5*(tempT.dot(temp)-I_matrix))
        self.strain = np.array(strain.copy())
        print(np.max(np.abs(self.strain[:,0,0])))
        print('strain calculated!!!')
        
    def clinicalStrainAxis(self, sampleCoord=None, order=2):
        '''
        clinical strain axis (direction): longitudinal strain, circumferential strain, radial strain, area strain
        '''
        # need test
        if type(sampleCoord)==type(None):
            sampleCoord=np.array(self.sampleCoord.copy())
        else:
            self.sampleCoord=np.array(sampleCoord.copy())
        innerFaceCenter = np.array(self.innerFaceCenter.copy())
        innerFaceNormal = np.array(self.innerFaceNormal.copy())
        outerFaceCenter = np.array(self.outerFaceCenter.copy())
        outerFaceNormal = np.array(self.outerFaceNormal.copy())
        
        innerDistValue=[]
        innerDistIndex=[]
        outerDistValue=[]
        outerDistIndex=[]
        sampleCoord = sampleCoord[sampleCoord[:,self.longAxis]>self.longAxisMin]
        sampleCoord = sampleCoord[sampleCoord[:,self.longAxis]<self.longAxisMax]
        self.sampleCoord = sampleCoord.copy()
        for i in sampleCoord:
            sample=i.reshape((1,3))
            #sample=np.append(sample,sample,axis=0)
            innerValue,innerIndex = np.array(self.distCalc(innerFaceCenter,sample,mode='min'))
            outerValue,outerIndex = np.array(self.distCalc(outerFaceCenter,sample,mode='min'))
            innerDistValue.append(innerValue)
            innerDistIndex.append(innerIndex)
            outerDistValue.append(outerValue)
            outerDistIndex.append(outerIndex)
        # if memory big enough, use codes below
        #innerDistValue,innerDistIndex = np.array(self.distCalc(innerFaceCenter,sampleCoord,mode='min'))
        #outerDistValue,outerDistIndex = np.array(self.distCalc(outerFaceCenter,sampleCoord,mode='min'))
        innerDistIndex = np.array(innerDistIndex).astype(int).reshape((1,-1),order='C')
        innerDistIndex = innerDistIndex.tolist()[0]
        innerDistValue = np.array(innerDistValue).reshape((1,-1),order='C').transpose()
        outerDistIndex = np.array(outerDistIndex).astype(int).reshape((1,-1),order='C')
        outerDistIndex = outerDistIndex.tolist()[0]
        outerDistValue = np.array(outerDistValue).reshape((1,-1),order='C').transpose()
        
        innerSampleNormal = innerFaceNormal[innerDistIndex]
        outerSampleNormal = outerFaceNormal[outerDistIndex]
        innerWeight=outerDistValue**order/(innerDistValue**order+outerDistValue**order)
        outerWeight=innerDistValue**order/(innerDistValue**order+outerDistValue**order)
        sampleNormal=(-1)*innerSampleNormal*innerWeight+outerSampleNormal*outerWeight
        normalize_v3(sampleNormal)
        self.sampleNormal=sampleNormal.copy()
        
        centerLine = self.centerLine.copy()
        self.radialAxis = sampleNormal.copy()
        n = np.cross( centerLine[1,:] - centerLine[0,:]  , self.radialAxis )
        normalize_v3(n)
        self.circumAxis = n.copy()
        n = np.cross( self.radialAxis  , self.circumAxis )
        normalize_v3(n)
        self.longitAxis = n.copy()
        
    def clinicalStrainCalc(self):
        strain = self.strain.copy()
        shape = strain.shape
        longitAxis = self.longitAxis.copy().reshape((shape[0],1,shape[2]))
        circumAxis = self.circumAxis.copy().reshape((shape[0],1,shape[2]))
        radialAxis = self.radialAxis.copy().reshape((shape[0],1,shape[2]))
        
        temp0=[]
        temp1=[]
        temp2=[]
        eigenVector = []
        eigenValue = []
        for i in range(shape[0]):
            # need to improve the algorithm, not use for loop
            temp0.append(longitAxis[i].dot(strain[i].dot(longitAxis[i].transpose())))
            temp1.append(circumAxis[i].dot(strain[i].dot(circumAxis[i].transpose())))
            temp2.append(radialAxis[i].dot(strain[i].dot(radialAxis[i].transpose())))
            w,v = np.linalg.eig(strain[i])      #one eigenvector per column in v
            eigenValue.append(w)
            eigenVector.append(v)
        
        self.longitStrain = np.array(temp0).reshape((shape[0],1))
        self.circumStrain = np.array(temp1).reshape((shape[0],1))
        self.radialStrain = np.array(temp2).reshape((shape[0],1))
        self.eigenValue = np.array(eigenValue)
        self.eigenVector = np.array(eigenVector)
        
        print('function clinicalStrainCalc done ^_^')
