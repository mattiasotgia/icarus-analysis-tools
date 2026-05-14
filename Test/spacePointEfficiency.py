# Test file
import os,sys
import math
import numpy  as np
import pandas as pd
from scipy.optimize import curve_fit
import plotly.graph_objects as go
import plotly.subplots as subplots
from plotly.colors import DEFAULT_PLOTLY_COLORS

import ROOT

import galleryUtils

print("name:",__name__)

import SBNDservices as services

# switch to "old" SBND geometry
#services.ServiceManager.setConfiguration("simulationservices_altgeometry_sbnd.fcl","sbnd_simulation_services")

configNew,tableNew = services.ServiceManager.defaultConfiguration()

print("Config:",configNew,", table:",tableNew)

# Add the path to our functions
sigProcPath = "/home/usher/test/icarus-analysis-tools"
sys.path.insert(0,sigProcPath)

from analysis_functions.getServices import getServices

detClocks,detProperties,geometryCore,lar_properties = getServices(services)

from analysis_functions.badChannels import getSBNDBadChannels

badChannels = getSBNDBadChannels()

offsets = [0., 0., 0.]

from analysis_functions.relateObjects import *
from analysis_functions.analyzeObjects import *


# define a constant for converting from number of electrons to ADC counts
nElectronsToADC = 1./110.

inputFile = "/media/usher/Drive2/sbnddata/decode_data_evb01_EventBuilder1_art1_run15772_10_20240803T213534-8a27789a-700f-400e-b499-502608b17270_Reco1-20240815T193145.root"
inputFile = "/media/usher/Drive2/sbnddata/prodmpvmpr_sbnd_MPVMPR-20241116T173006_G4-20241116T173726_DetSim-20241116T175023_6c239230-2e28-444f-8561-46c411582c3f_Reco1-20241118T152246.root"
inputFile = "/media/usher/Drive2/sbnddata/prodmpvmpr_sbnd_MPVMPR-20241116T173006_G4-20241116T173726_DetSim-20241116T175023_6c239230-2e28-444f-8561-46c411582c3f_Reco1-20241118T152246_Reco1Test-20241204T182718.root"
inputFile = "/media/usher/Drive2/sbnddata/prodmpvmpr_sbnd_MPVMPR-20241116T173006_G4-20241116T173726_DetSim-20241116T175023_6c239230-2e28-444f-8561-46c411582c3f_Reco1-20241118T152246_Reco1Test-20241205T194018.root"
inputFile = "/media/usher/Drive2/sbnddata/prodmpvmpr_sbnd_MPVMPR-20241206T174559_G4-20241206T175616_DetSim-20241206T181328_8935f56f-ac6e-4939-a7d1-ae5e5e3cf363_Reco1-20241206T195725.root"
#inputFile = "/media/usher/Drive2/sbnddata/prodmpvmpr_sbnd_MPVMPR-20241206T174559_G4-20241206T175616_DetSim-20241206T181328_8935f56f-ac6e-4939-a7d1-ae5e5e3cf363_Reco1-20241206T195725_Reco1Test-20241211T190048.root"
#inputFile = "/media/usher/Drive2/sbnddata/prodmpvmpr_sbnd_MPVMPR-20241206T174559_G4-20241206T175616_DetSim-20241206T181328_8935f56f-ac6e-4939-a7d1-ae5e5e3cf363_Reco1-20241206T195725_Reco1Test-20241212T164217.root"
#inputFile = "/media/usher/Drive2/sbnddata/prodmpvmpr_sbnd_MPVMPR-20241116T173006_G4-20241116T173726_DetSim-20241116T175023_6c239230-2e28-444f-8561-46c411582c3f_Reco1-20241212T173428.root"
inputFile = "/media/usher/Drive2/sbnddata/fixed/prodmpvmpr_sbnd_MPVMPR-20241212T175410_G4-20241212T175846_DetSim-20241212T180653_cd83ad1c-9bc8-4059-8c22-f5fd3ce377df_Reco1-20241216T174835.root"
#inputFile = "/media/usher/Drive2/sbnddata/bear1ddecon/prodmpvmpr_sbnd_MPVMPR-20250114T174306_G4-20250114T175650_DetSim-20250114T175901_Reco1-20250114T181047.root"
#inputFile = "/media/usher/Drive2/sbnddata/bear1ddecon/prodmpvmpr_sbnd_MPVMPR-20250114T174306_G4-20250114T221255_DetSim-20250114T221535_Reco1-20250114T221637.root"
inputFile = "/media/usher/Drive2/sbnddata/v09_72_01/prodsingle_sbnd_SinglesGen-20250125T003706_G4-20250125T005357_DetSim-20250125T012007_Reco1-20250127T173357_supera-20250130T192202.root"
#inputFile = "/media/usher/Drive2/sbnddata/reco1_doublets_hitcfg.root"   # This has doublets but also bad channels + threshold over baseline for hit finding
#inputFile = "/media/usher/Drive2/sbnddata/intime/prodmpvmpr_sbnd_MPVMPR-20250131T184911_G4-20250131T185426_DetSim-20250201T000632_Reco1-20250201T001017.root" # As above but intime
#inputFile = "/media/usher/Drive2/sbnddata/particlegun/prodsingle_sbnd_SinglesGen-20250211T194135_G4-20250211T195454_DetSim-20250211T210923_Reco1-20250211T220356.root" # As above but intime
inputFile = "/media/usher/Drive2/sbnddata/particlegun/prodsingle_sbnd_SinglesGen-20250211T194135_G4-20250211T195454_DetSim-20250212T210702_Reco1-20250212T211432.root" # As above but intime
inputFile = "/media/usher/Drive2/sbnddata/reco1-art.root" # Francois bad event
#inputFile = "/media/usher/Drive2/sbnddata/reco1-art_Redo-20250226T183846.root" # above but re-run cluster3d

#sampleEvents = galleryUtils.makeEvent(inputFile)
#sampleEvents = gallery::Event(inputFile)
#sampleEvents = ROOT.gallery.Event(galleryUtils.makeFileList([ ROOT.std.string(inputFile) ]))

#####################################################################################################################
#
# Loop over events
#

inputFiles = ["/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175410_G4-20241212T175846_DetSim-20241212T180653_cd83ad1c-9bc8-4059-8c22-f5fd3ce377df_Reco1-20241216T174835.root",
              "/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175350_G4-20241212T175854_DetSim-20241212T180731_ae089d32-8d6e-4c2f-9f59-63a8705f26ec_Reco1-20241216T175024.root",
              "/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175426_G4-20241212T180425_DetSim-20241212T181748_1bc74075-9b78-49cc-b3fc-65b3f2894810_Reco1-20241216T181918.root",
              "/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175416_G4-20241212T180403_DetSim-20241212T181701_4c53bb70-9b9a-4671-8c47-735949a27043_Reco1-20241216T182118.root",
              "/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175357_G4-20241212T175937_DetSim-20241212T180806_82183a4f-3cfa-4b8f-bc14-ecedd1e51749_Reco1-20241216T182258.root",
              "/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175422_G4-20241212T175911_DetSim-20241212T180846_fe62972e-dd71-403b-8190-729e60e2cb40_Reco1-20241216T183118.root",
              "/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175426_G4-20241212T175934_DetSim-20241212T180841_7e53f8c7-907a-4f0c-aa09-f223ff1f48d7_Reco1-20241216T183542.root",
              "/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175445_G4-20241212T180125_DetSim-20241212T181054_e7464944-ac08-458e-919d-e556a4cd0097_Reco1-20241216T183715.root", 
              "/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175445_G4-20241212T180117_DetSim-20241212T181038_929f3944-6149-4d0d-8f19-f2ee90240620_Reco1-20241216T183838.root", 
              "/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175445_G4-20241212T180140_DetSim-20241212T181112_b9620255-9bc8-4021-9dc9-9ba20aaddc0b_Reco1-20241216T184048.root",
              "/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175426_G4-20241212T175927_DetSim-20241212T180839_c2a4b020-6a7c-4380-bde9-cd55f2f6ad4b_Reco1-20241216T184514.root",
              "/media/usher/Drive2/sbnddata/files1216/prodmpvmpr_sbnd_MPVMPR-20241212T175445_G4-20241212T180128_DetSim-20241212T181102_60f04581-ac97-4000-8627-ffbdf7fe3ce7_Reco1-20241216T184958.root" ]

# If just using a single input file try this
inputFiles   = galleryUtils.makeFileList(inputFile)

# Loop over events in the input file looking for the one we want
#for event in galleryUtils.forEach(sampleEvents):
spDataTag      = "cluster3d" #::Reco1Test"
hitAssnProdTag = "cluster3d" #::Reco1Test"
hitDataTag     = "cluster3d"
hitDataTag     = "gaushit"
simChannelTag  = "simtpc2d:simpleSC"
#simChannelTag  = "simdrift"
mcParticleTag  = "simplemerge"
simEnergyTag   = "largeant:LArG4DetectorServicevolTPCActive"
sedLiteTag     = "sedlite"

# Register with gallery
sampleEvents = ROOT.gallery.Event(inputFiles)

for iEvent, event in enumerate(galleryUtils.forEach(sampleEvents)):    
    # Recover handles to our data objects
    getHitHandle        = sampleEvents.getValidHandle[ROOT.vector[ROOT.recob.Hit]](hitDataTag)
    getSpacePointHandle = sampleEvents.getValidHandle[ROOT.vector[ROOT.recob.SpacePoint]](spDataTag)
    getSimChannelHandle = sampleEvents.getValidHandle[ROOT.vector[ROOT.sim.SimChannel]](simChannelTag)
    #getSimEnergyHandle  = sampleEvents.getValidHandle[ROOT.vector[ROOT.sim.SimEnergyDeposit]](simEnergyTag)
    
    # Convert these to the actual collections
    simChannelCol       = getSimChannelHandle.product()
    gausHitCol          = getHitHandle.product()
    spacePointCol       = getSpacePointHandle.product()
                
    # Associations
    hitAssns = sampleEvents.getValidHandle[ROOT.art.Assns[ROOT.recob.Hit,ROOT.recob.SpacePoint]](hitAssnProdTag).product()

    print("simChannelCol len:",len(simChannelCol),", gausHitCol len:",len(gausHitCol),", spacePointCol len:",len(spacePointCol))
    print("Size of hit<->SpacePoint assns:",hitAssns.size())

    # Now recover track/space point/hit info
    trackToIDEMap      = getTrackToIDEMap(simChannelCol,detClocks,badChannels)
    #trackToMCPartMap   = getTrackToMCPartMap(mcParticleCol)
    channelToHitMap    = getChannelToHitMap(gausHitCol)
    #channelToSPHitMap  = getChannelToHitMap(spHitCol)

    print("Return from channelToHitMap:",len(channelToHitMap),", tracktoIDEMap has:",len(trackToIDEMap))
    print("Calling getTrackToHitMap")
    
    trackToHitMap,trackToHitPeakIdeMap,trackToNoChanMap,trackToMissHitMap = getTrackToHitMap(trackToIDEMap,channelToHitMap)

    print("Returned from getTrackToHitMap, len:",len(trackToHitMap))
    
    #spPtrDict,hitPtrList,spacePointDict,hitToSPDict,positionDict,hitSPList = getSpacePointInfo(hitAssns)
    spPtrDict,hitPtrList,spacePointDict,hitToSPDict,hitSPList = getSpacePointInfo(hitAssns)
    
    print("Returned from getSpacePointIno, len(trackToHitMap):",len(trackToHitMap),", hitToSPDict:",len(hitToSPDict))

    # Recover the sim energy deposit positions
    simEnergyCol = sampleEvents.getValidHandle[ROOT.vector[ROOT.sim.SimEnergyDeposit]](simEnergyTag).product()
    simPositions = {}
    negPositions = {}

    print("Length of SimEnergyDeposits:",len(simEnergyCol))

    sedCounter = 0
    
    for simIdx,simEnergyDeposit in enumerate(simEnergyCol) :
        trackID = simEnergyDeposit.TrackID()

        if sedCounter < 100:
            print("SED idx:",simIdx,", trackID:",trackID,", has time:",simEnergyDeposit.Time(),", step length:",simEnergyDeposit.StepLength())

        sedCounter += 1
       
        if trackID > 0:
            simPositions.setdefault(trackID,[]).append([simEnergyDeposit.X(),simEnergyDeposit.Y(),simEnergyDeposit.Z()])
        else:
            negPositions.setdefault(trackID,[]).append([simEnergyDeposit.X(),simEnergyDeposit.Y(),simEnergyDeposit.Z()])
        
    print("Sim Energy Deposits:",len(simPositions))

    for trackID,positionList in simPositions.items():
        if len(positionList) > 2:
            print("TrackID:",trackID,", len simPositions:",len(positionList))

    # Build the track to SpacePoint to Hit mapping
    trackToSPHitListMap = trackToSpacePointMap(trackToHitMap,hitToSPDict)
                    
    print("Len of trackToHitMap:",len(trackToHitMap),", and for trackToSPHitListMap:",len(trackToSPHitListMap))

    # These are the usual ipython objects, including this one you are creating
    ipython_vars = ['In', 'Out', 'exit', 'quit', 'get_ipython', 'ipython_vars']

    # Get a sorted list of the objects and their sizes
    #ilist = sorted([(x, sys.getsizeof(globals().get(x))) for x in dir() if not x.startswith('_') and x not in sys.modules and x not in ipython_vars], key=lambda x: x[1], reverse=True)
    ilist = sorted([(x, sys.getsizeof(globals().get(x))) for x in dir()], key=lambda x: x[1], reverse=True)
    print(ilist)

    totalList = sum(x[1] for x in ilist)
    print("Total list:",totalList)
    
    # Now let's see if we can associate space points that are associated to a track with 
    for trackID,spHitListMap in trackToSPHitListMap.items():
        print("--> Looping over tracks, trackID:",trackID,", len sphit map:",len(spHitListMap)) #,", simEnergyDeposits:",len(simPositions[trackID]))
        
        if trackID in simPositions:
            print("    with length of simPositions:",len(simPositions[trackID]))
        
            simEnergyPos = np.array(simPositions[trackID])

            print("    and simEnergyPos shape:",simEnergyPos.shape)
            
            countEm = 0
            
            for spacePointIdx,hitList in spHitListMap.items():
                if spacePointIdx not in spPtrDict:
                    print("Space point idx:",spacePointIdx," not in spPtrDict")
                    continue

                spacePoint    = spPtrDict[spacePointIdx]
                spacePointPos = np.array([spacePoint.position().X(),spacePoint.position().Y(),spacePoint.position().Z()])
                
                distances = np.sqrt(np.sum((simEnergyPos - spacePointPos)**2,axis=1))
                
                nearestIdx = np.argmin(distances)
                
                print("SP Pos:",spacePointPos,", nearest point:",simEnergyPos[nearestIdx],", distance:",distances[nearestIdx])
                
                if countEm > 20:
                    break
                    
                countEm += 1
    break
        

