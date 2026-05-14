# Note that the analysis functions below need to decode art root data objects so we will need 
# access to the gallery utilities to access theb
# All of this assumes you have set up the icarusalg package so you can interface to LArSoft services
#import ROOT
#import galleryUtils
import numpy as np

# Define a gaussian function
# x are the input coordinates
# PeakAmp is the amplitude of the gaussian
# PeakTime is the center of the gaussian
# sigma is the sigma of the gaussian
# baseline is the default baseline of the waveform
def gaus(x,peakAmp,peakTime,sigma,baseline):
       return peakAmp*np.exp(-0.5*(x-peakTime)**2/(sigma**2))+baseline


# Define a generic function to map channels to the individual objects in that channel
def mapChannelsToObjects(objectTags, getObjectHandle):
    # First task is to build the channel to ROI map
    channelToObjects = {}
    
    for objectTag in objectTags:
        objectCol = getObjectHandle(objectTag).product()
        
        for object in objectCol :
            channel = object.Channel()
            
            if channel in channelToObjects:
                channelToObjects[channel].append(object)
            else:
                channelToObjects[channel] = [object]
    
    return channelToObjects

# define a constant for converting from number of electrons to ADC counts
nElectronsToADC = 1./110.

# Provide a mapping between channel and hits on that channel
def getChannelToHitMap(hitCol):    
    channelToHitMap = {}

    for hit in hitCol :
        channelToHitMap.setdefault(hit.Channel(),[]).append(hit)
            
    return channelToHitMap

# Provide a mapping from MC TrackID to associated SimChannel IDEs
def getTrackToIDEMap(simChannelCol,detClocks,badChannels=[]):
    trackToIDEMap = {}

    for simChannel in simChannelCol:
        channel = simChannel.Channel()
        
        # skip bad channels
        if channel in badChannels:
            continue
            
        
        tdcides = simChannel.TDCIDEMap()
        
        for tdcide in tdcides:
            tick = detClocks.DataForJob().TPCTDC2Tick(tdcide.first)
            
            for ide in tdcide.second:
                if ide.trackID not in trackToIDEMap:
                    trackToIDEMap[ide.trackID] = {}
                    
                channelToIDEMap = trackToIDEMap[ide.trackID]
                
                if channel not in channelToIDEMap:
                    channelToIDEMap[channel] = {}
                    
                tickToIDEMap = channelToIDEMap[channel]
                
                tickToIDEMap.setdefault(int(tick),[]).append(ide)

    return trackToIDEMap

# Provide a mapping between MC TrackIDs and MCParticles
def getTrackToMCPartMap(mcParticleCol):
    trackToMCPartMap = {}

    for mcParticle in mcParticleCol:
        trackToMCPartMap[mcParticle.TrackId()] = mcParticle

    return trackToMCPartMap
    
# Mega fun ction to build space point to/from hit associations 
def getSpacePointInfo(hitAssns):
    spPtrDict      = {}
    hitPtrList     = []
    spacePointDict = {}
    hitToSPDict    = {}
    #hitToPosDict   = {}
    hitSPList      = []

    # Build up dictionaires between space points and hits, etc. 
    for hitPtr,spacePointPtr in hitAssns:
        if hitPtr is None or spacePointPtr is None:
            print("Null pointers for hit assocations")
            continue

        #spacePointDict.setdefault((spacePointPtr.id(),spacePointPtr.key()),[]).append(hitPtr)
        spPtrDict[spacePointPtr.ID()] = spacePointPtr.get()  # There will be no duplicates in the dictionary
        hitPtrList.append(hitPtr.get())                      # Can be duplicates here 
        spacePointDict.setdefault(spacePointPtr.ID(),[]).append(hitPtr.get())
        hitToSPDict.setdefault(hitPtr.get(),[]).append(spacePointPtr.ID())
        #hitToPosDict[hitPtr.get()] = spacePointPtr.position()
        
    # Trim out duplicates
    hitPtrList = list(set(hitPtrList))
    
    # Go thgrough again to pick out hits in 2 or 3 hit spacepoints
    for key,hitList in spacePointDict.items():
        if len(hitList) < 3:   # change this from!= 3:
            continue
            
        for hit in hitList:
            hitSPList.append(hit)
            
    # Remove duplicates
    hitSPList = list(set(hitSPList))
    
    #return spPtrDict,hitPtrList,spacePointDict,hitToSPDict,hitToPosDict,hitSPList
    return spPtrDict,hitPtrList,spacePointDict,hitToSPDict,hitSPList

# Basically we are building a map from MC Track ID to associated hits
# We also need to build in a "miss" map too...
def getTrackToHitMap(trackToIDEMap,channelToHitMap):
    trackToHitMap        = {}
    trackToNoChanMap     = {}
    trackToMissHitMap    = {}
    trackToHitPeakIdeMap = {}
    
    # Loop through the trackIDs to start
    for trackID,channelToTickIDEs in trackToIDEMap.items():
        # Now through each channel with ides associated to this MC Track
        trackToHitPeakIdeMap.setdefault(trackID,{})
        hitToPeakIdeMap = trackToHitPeakIdeMap[trackID]
        
        for channel,tickToIDEMap in channelToTickIDEs.items():
            #Set up to find the tick with the peak number of electrons that is associated to this track
            peakTick = 0
            peakElec = 0
            peakIde  = None
            
            for tick,ides in tickToIDEMap.items():
                peakElectrons = 0
                trackIde = None
                for ide in ides:
                    if ide.trackID == trackID:
                        peakElectrons = ide.numElectrons
                        trackIde = ide
                        break
                        
                if peakElectrons > peakElec:
                    peakTick = tick
                    peakElec = peakElectrons
                    peakIde  = trackIde
            
            peakElec *= nElectronsToADC
                
            # Put a lower bound on the peak electrons (in ADC equivalents)
            if peakElec < 8:
                continue
            
            # Now find the "best" associated hit (best defined as closest)
            if channel in channelToHitMap:
                hitList = channelToHitMap[channel]

                hitList.sort(key=lambda hit: hit.PeakTime())
                
                localHitList = []
                
                # Loop over hits, note special handling for "long hits"
                for hitIdx in range(len(hitList)):
                    hit = hitList[hitIdx]

                    #print("-- Channel:",channel,", hitIdx:",hitIdx,", dog:",hit.DegreesOfFreedom(),", index:",hit.LocalIndex())
                    
                    if hit.DegreesOfFreedom() > 1.:
                        localHitList.append([hit])
                    else:
                        if hit.LocalIndex() == 0 or len(localHitList) < 1:
                            localHitList.append([hit])
                        else:
                            localHitList[-1].append(hit)
               
                bestDeltaPeakTime = -10000
                peakAmp           = 0
                bestHitIdx        = -1
                
                # Loop again handling the separate types of hits
                for hitsIdx,hits in enumerate(localHitList):
                    if len(hits) > 1:
                        startTick = hits[0].StartTick()
                        endTick   = hits[0].EndTick()
                        peakTime  = startTick + (endTick - startTick)/2
                        peakSigma = (endTick - startTick) / 5.
                    else:
                        peakTime  = hits[0].PeakTime()
                        peakSigma = hits[0].RMS()
                        
                    deltaPeak = peakTime - peakTick
                        
                    if abs(deltaPeak) < 5.*peakSigma and abs(deltaPeak) < abs(bestDeltaPeakTime):
                        #if abs(deltaPeak) < abs(bestDeltaPeakTime):
                        bestDeltaPeakTime = deltaPeak
                        bestHitIdx        = hitsIdx
                        
                if bestHitIdx > -1:
                    for hit in localHitList[bestHitIdx]:                    
                        trackToHitMap.setdefault(trackID,[]).append(hit)
                        hitToPeakIdeMap[hit] = [peakTick,peakIde]
                else:
                    trackToMissHitMap.setdefault(trackID,[]).append([channel,peakTick,peakIde])
            else:
                trackToNoChanMap.setdefault(trackID,[]).append(channel)
                    
    return trackToHitMap,trackToHitPeakIdeMap,trackToNoChanMap,trackToMissHitMap

# Here we build the mapping between MC TrackID and space points
def trackToSpacePointMap(trackToHitMap,hitToSPDict):
    # The strategy is to loop through hits associated to trackIDs and build a frequency mapping to space points
    trackToSPHitListMap = {}

    trackSPHitList = hitToSPDict.keys()
    
    for trackID,hitList in trackToHitMap.items():
        if trackID not in trackToSPHitListMap:
            trackToSPHitListMap[trackID] = {}
            
        spHitListMap = trackToSPHitListMap[trackID]

        countem = 0
        
        for hit in hitList:
            # Note that the "==" operator is not defined for hits so need to compare channel # and PeakTime
            matchedHit = next((x for x in trackSPHitList if x.Channel() == hit.Channel() and abs(x.PeakTime()-hit.PeakTime()) < 0.01), None)

            if trackID == -2:
                #if matchedHit is not None:
                    #print("  - count",countem,", hit channel/time:",hit.Channel(),"/",hit.PeakTime(),", matchedHit:",matchedHit.Channel(),"/",matchedHit.PeakTime())
                #else:
                if matchedHit is None:
                    print("  - count",countem,", matchedHit not found for hit channel/time:",hit.Channel(),"/",hit.PeakTime())
            countem += 1
            
            if matchedHit is not None:
                for spacePointIdx in hitToSPDict[matchedHit]:
                    spHitListMap.setdefault(spacePointIdx,[]).append(matchedHit)
    
    return trackToSPHitListMap
