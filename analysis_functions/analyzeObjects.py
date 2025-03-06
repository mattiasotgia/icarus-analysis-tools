# Note that the analysis functions below need to decode art root data objects so we will need 
# access to the gallery utilities to access theb
# All of this assumes you have set up the icarusalg package so you can interface to LArSoft services
import numpy as np


def computePCA(inputPoints,iterate=False,nSigma=3):
    # make a copy
    pointCloud = inputPoints.copy()
    
    # Compute PCA
    #print("computePCA, pointCloud shape:",pointCloud.shape)
    cloudMean   = pointCloud.mean(axis=0)
    localCloud  = pointCloud - cloudMean
    #print("computePCA, cloudMean:",cloudMean)
    #pointCov    = pointCloud.T @ pointCloud / xRange.size
    pointCov    = np.cov((localCloud).T)
    
    eigenVals   = np.zeros(3)
    eigenVecs   = np.zeros((3,3))
    
    if np.any(np.isinf(pointCov)) or np.any(np.isnan(pointCov)):
        print("Covariance matrix has bad values:")
    else:
        eigenVals,eigenVecs = np.linalg.eig(pointCov)
        
        # Iterate to throw out outliers?
        #if iterate:
        #    if eigenVals[1] > eigenVals[0]:
        #        eigenVals = np.flip(eigenVals)
        #        eigenVecs = np.flip(eigenVecs,axis=1)
        #                
        #    if eigenVecs[0,0] < 0.:
        #        eigenVecs *= -1.
        #        
        #    eigenVecs[:,1] = np.cross(eigenVecs[:,0],eigenVecs[:,1]) * eigenVecs[:,1]
        #    
        #    rangeToKeep   = nSigma * np.sqrt(eigenVals[1])
        #    rotPointCloud = (eigenVecs.T @ localCloud.T).T
        #    points        = pointCloud[np.where(np.absolute(rotPointCloud[:,1])<rangeToKeep)]

        #    if np.size(points,axis=0) < np.size(pointCloud,axis=0) and np.size(points,axis=0) > 6:
        #        eigenVals,eigenVecs,points,cloudMean = computePCA(points[:,0],points[:,1],iterate,nSigma)
            
        pointCloud -= cloudMean
        
    return eigenVals,eigenVecs,pointCloud,cloudMean
