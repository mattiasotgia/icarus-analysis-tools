import numpy as np

import plotly.graph_objects as go
import plotly.subplots as subplots
from   plotly.colors import DEFAULT_PLOTLY_COLORS

# Useful for plotting...
def gaus(x,peakAmp,peakTime,sigma,baseline):
    return peakAmp*np.exp(-0.5*(x-peakTime)**2/(sigma**2))+baseline
    
def box(x,peakAmp,peakTime,sigma,baseline):
    boxArray     = np.full(len(x),peakAmp) + baseline
    boxArray[0]  = baseline
    boxArray[-1] = baseline
    return boxArray
    
def plotWireROIs2D(wireROICol, wireIDs, geometryCore, offsets, colors, fig):
    print("Plotting Wire data...")
    roiTimeVals = []
    roiWireVals = []
    
    for wireROI in wireROICol :
        wireIDVec = geometryCore.ChannelToWire(wireROI.Channel())
        wireIDidx = -1
        if wireIDVec[0] == wireIDs[0]:
            wireIDidx = 0
        elif wireIDVec.size() > 1 and wireIDVec[1] == wireIDs[0]:
            wireIDidx = 1
            
        if wireIDidx > -1:
            #print("Channel:",channelROI.Channel(),", size:",channelROI.NSignal())
            dataVals  = np.zeros(4096)
            for roi in wireROI.SignalROI().get_ranges() :
                startTick = roi.begin_index() - offsets[wireIDVec[wireIDidx].Plane]
                maxTick   = np.argmax(np.array(roi.data())) + startTick
                stopTick  = startTick + roi.data().size()
                wireLow   = wireIDVec[wireIDidx].Wire - 0.5
                fig.add_shape(type="rect",x0=wireLow, y0=startTick, x1=wireLow+1., y1=stopTick,line=dict(color=colors))
                
                roiTimeVals.append(maxTick)
                roiWireVals.append(wireIDVec[wireIDidx].Wire)

    fig.add_trace(go.Scatter(x=roiWireVals, y=roiTimeVals, mode='markers'))
    fig.update_traces(marker=dict(size=4,line=dict(width=1,color=colors)))
          
    print("Done plotting WireROI data")
                                   
    return
    
def plotChannelROIs2D(channelROICol, wireIDs, geometryCore, offsets, color, fig):
    print("Plotting Wire data...")
    roiTimeVals = []
    roiWireVals = []

    rectsToAdd  = []
    
    for channelROI in channelROICol :
        wireIDVec = geometryCore.ChannelToWire(channelROI.Channel())
        wireIDidx = -1
        if wireIDVec[0].asPlaneID() == wireIDs[0].asPlaneID():
            wireIDidx = 0
        elif wireIDVec.size() > 1 and wireIDVec[1].asPlaneID() == wireIDs[0].asPlaneID():
            wireIDidx = 1
            
        if wireIDidx > -1:
            for roi in channelROI.SignalROIF().get_ranges() :
                startTick = roi.begin_index() - offsets[wireIDVec[wireIDidx].Plane]
                maxTick   = np.argmax(np.array(roi.data())) + startTick
                stopTick  = startTick + roi.data().size() + 1
                wireLow   = wireIDVec[wireIDidx].Wire - 0.5
                rectsToAdd.append(dict(type="rect",x0=wireLow, y0=startTick, x1=wireLow+1., y1=stopTick,line=dict(color=color)))
                #fig.add_shape(type="rect",x0=wireLow, y0=startTick, x1=wireLow+1., y1=stopTick,line=dict(color=color))
                
                roiTimeVals.append(maxTick)
                roiWireVals.append(wireIDVec[wireIDidx].Wire)

    fig.add_trace(go.Scatter(x=roiWireVals, y=roiTimeVals, mode='markers',name="ROI Max Ticks"))
    fig.update_traces(marker=dict(size=4,line=dict(width=1,color=color)))
    fig.update_layout(shapes=rectsToAdd)
    fig.update_shapes(name="ROI regions")
          
    print("Done plotting ChannelROI data")
                                   
    return
    
def plotWireROIs(wireROICol, wireIDs, geometryCore, offsets, colors, fig):
    print("Plotting Wire data...")
    dataVals = np.zeros(4096)
    
    for wireROI in wireROICol :
        wireIDVec = geometryCore.ChannelToWire(wireROI.Channel())
        wireIDidx = -1
        if wireIDVec[0] == wireIDs[0]:
            wireIDidx = 0
        elif wireIDVec.size() > 1 and wireIDVec[1] == wireIDs[0]:
            wireIDidx = 1
            
        if wireIDidx > -1:
            #print("Channel:",wireROI.Channel(),", size:",wireROI.NSignal())

            for roi in wireROI.SignalROI().get_ranges() :
                startTick = roi.begin_index() - offsets[wireIDVec[wireIDidx].Plane]
                stopTick  = startTick + roi.data().size()

                dataVals[startTick:stopTick] = roi.data()

    label   = str('WireROI Channel: %d' % (wireROI.Channel()))
    fig.add_trace(go.Scatter(x=np.arange(0,4096), y=dataVals,name=label,line_shape="hvh",line=dict(color=colors,width=1)),row=wireIDs[0].Plane+1,col=1)
    #fig.update_traces(marker=dict(size=4,line=dict(width=1,color=colors)),row=wireIDs[0].Plane+1,col=1)
          
    print("Done plotting WireROI data")
                                   
    return
    
def plotChannelROIs(channelROICol, wireIDs, geometryCore, offsets, colors, fig):
    print("Plotting Wire data...")
    dataVals = np.zeros(4096)
    
    for channelROI in channelROICol :
        wireIDVec = geometryCore.ChannelToWire(channelROI.Channel())
        wireIDidx = -1
        if wireIDVec[0] == wireIDs[0]:
            wireIDidx = 0
        elif wireIDVec.size() > 1 and wireIDVec[1] == wireIDs[0]:
            wireIDidx = 1
            
        if wireIDidx > -1:
            for roi in channelROI.SignalROIF().get_ranges() :
                startTick = roi.begin_index() - offsets[wireIDVec[wireIDidx].Plane]
                stopTick  = startTick + roi.data().size()

                dataVals[startTick:stopTick] = roi.data()

    label   = str('ChannelROI Channel: %d' % (channelROI.Channel()))
    fig.add_trace(go.Scatter(x=np.arange(0,4096), y=dataVals,name=label,line_shape="hvh",line=dict(color=colors[wireIDs[0].Plane],width=1)),row=wireIDs[0].Plane+1,col=1)
    #fig.update_traces(marker=dict(size=4,line=dict(width=1,color=colors)),row=wireIDs[0].Plane+1,col=1)
          
    print("Done plotting ChannelROI data")
                                   
    return
    
def plotRawDigits(rawDigitCol, wireIDs, geometryCore, colors, fig):
    print("Plotting RawDigit data...")
    
    for wireID in wireIDs:
        channel = geometryCore.PlaneWireToChannel(wireID)
        
        print("plotRawDigits looking for channel:",channel," (from TPC/Plane/Wire:",wireID.TPC,"/",wireID.Plane,"/",wireID.Wire)
        
        for rawDigit in rawDigitCol:
            if rawDigit.Channel() == channel:
                rawDigitVals = np.array(rawDigit.ADCs()) - rawDigit.GetPedestal()
                label   = str('RawDigit Channel: %d' % (channel))
                fig.add_trace(go.Scatter(x=np.arange(0,4096),y=rawDigitVals,name=label,line_shape="hvh",line=dict(color='black',width=1)),row=wireID.Plane+1,col=1)
                print("Found channel with pedestal:",rawDigit.GetPedestal())
                
    print("Done plotting RawDigits")
                                   
    return
    
# Plot hits on a waveform
def plotHits(hitCol, wireIDs, geometryCore, color, label, lineStyle, fig):
    dataVals   = np.zeros(4096)
    xTicksList = []
    
    print("PlotHits looping over",len(hitCol)," hits, channels:",wireIDs)

    for hit in hitCol :
        wireIDVec = geometryCore.ChannelToWire(hit.Channel())
        wireIdx   = -1
        
        if wireIDVec[0].TPC == wireIDs[0].TPC and wireIDVec[0].Plane == wireIDs[0].Plane and wireIDVec[0].Wire == wireIDs[0].Wire:
            wireIdx = 0
        elif wireIDVec.size() > 1 and wireIDVec[1].TPC == wireIDs[0].TPC and wireIDVec[1].Plane == wireIDs[0].Plane and wireIDVec[1].Wire == wireIDs[0].Wire:
            wireIdx = 1
            
        if wireIdx > -1:
            wireID = wireIDVec[wireIdx]

            rmsFactor = 1.

            if hit.DegreesOfFreedom() > 1:            
                xTicks = np.arange(max(round(hit.PeakTimeMinusRMS(4.)-1.),0),min(round(hit.PeakTimePlusRMS(4.)+1.),4095))
                #print("xticks:",xTicks)
                hitVals = gaus(xTicks,hit.PeakAmplitude(),round(hit.PeakTime()),hit.RMS(),0)
    
                xMinVal = int(xTicks[0])
                xMaxVal = xMinVal+len(xTicks)
                dataVals[xMinVal:xMaxVal] += hitVals
                print("-Channel:",hit.WireID(),", min/max:",xMinVal,"/",xMaxVal)
                xTicksList.append(xTicks)
            else:
                peakTime = hit.PeakTime()
                peakRMS  = hit.RMS() #/2.
                xTicks = np.arange(max(round(peakTime-3.*peakRMS),0),min(round(peakTime+3.*peakRMS),4095))
                #print("xticks:",xTicks)
                hitVals = box(xTicks,0.5*hit.PeakAmplitude(),round(hit.PeakTime()),hit.RMS(),0)
                
                xMinVal = int(xTicks[0])
                xMaxVal = xMinVal+len(xTicks)
                dataVals[xMinVal:xMaxVal] += hitVals
                print("-Channel:",hit.WireID(),", min/max:",xMinVal,"/",xMaxVal)
                xTicksList.append(xTicks)

                rmsFactor = 0.5

            # Plot the peak time and RMS
            fig.add_trace(go.Scatter(x=[hit.PeakTime(),hit.PeakTime()],y=[0,rmsFactor*hit.PeakAmplitude()],name="Hit Time: "+str(hit.PeakTime()),mode='lines',line=dict(color="gray",width=1)),row=wireIDs[0].Plane+1,col=1)
            fig.add_trace(go.Scatter(x=[hit.PeakTimeMinusRMS(1.),hit.PeakTimePlusRMS(1.)],y=[rmsFactor*hit.PeakAmplitude()/2.,rmsFactor*hit.PeakAmplitude()/2.],name="Hit RMS: "+str(hit.RMS()),mode='lines',line=dict(color="gray",width=1)),row=wireIDs[0].Plane+1,col=1)

     
    label = label+str(' Hits: %d' % wireIDs[0].Wire)
    
    lastX = -1
    for xTicks in xTicksList:
        startIdx  = 0
        xTicksLen = len(xTicks)
        for x in xTicks:
            if x <= lastX:
                startIdx  += 1
                xTicksLen -= 1
        if xTicksLen > 0:
            fig.add_trace(go.Scatter(x=xTicks[startIdx:],y=dataVals[int(xTicks[startIdx]):int(xTicks[startIdx])+xTicksLen],name=label,mode='lines',line=dict(color=color,width=1,dash=lineStyle)),row=wireIDs[0].Plane+1,col=1)
            
            
    return

# Plot SimChannel info as a waveform
def plotSimChannels(simChannelCol, wireIDs, wireGeom, detClocks, colors, fig):
    print("plotting sim for TPC:",wireIDs[0].TPC,", Plane:",wireIDs[0].Plane,", Wire:",wireIDs[0].Wire)

    for simChannel in simChannelCol:
        wireIDVec = wireGeom.ChannelToWire(simChannel.Channel())

        # There will only be 1 cryostat even if two WireIDs...
        if wireIDVec[0].Cryostat != wireIDs[0].Cryostat:
            continue
        
        wireIDidx = -1
        if wireIDVec[0].TPC == wireIDs[0].TPC and wireIDVec[0].Plane == wireIDs[0].Plane and wireIDVec[0].Wire == wireIDs[0].Wire:
            wireIDidx = 0
        elif wireIDVec.size() > 1 and wireIDVec[1].TPC == wireIDs[0].TPC and wireIDVec[1].Plane == wireIDs[0].Plane and wireIDVec[1].Wire == wireIDs[0].Wire:
            wireIDidx = 1
            
        if wireIDidx > -1:
            #print("Found SimChannel info:",simChannel)
            # loop through all ticks and recover charge
            wireID = wireIDVec[wireIDidx]
            print("SimChannel for",simChannel.Channel(),", Cryostat:",wireID.Cryostat,", TPC:",wireID.TPC,", Plane:",wireID.Plane,", Wire:",wireID.Wire,", wireIDidx:",wireIDidx)
            
            tdcides  = simChannel.TDCIDEMap()
            ideVals = np.zeros(4096)
            trackDict = {} #defaultdict(list)
            
            for tdcide in tdcides:
                tick = detClocks.DataForJob().TPCTDC2Tick(tdcide.first)

                if tick < 0 or tick > 4095:
                    continue
                #print("TDC value:",tdcide.first,", tick:",tick)

                for ide in tdcide.second:
                    #if ide.trackID < 0:
                    #    continue
                    if ide.trackID not in trackDict:
                        trackDict[ide.trackID] = np.zeros(4096)
                    #ideVals[int(tick)] += ide.numElectrons / 60.
                    trackDict[ide.trackID][int(tick)] += ide.numElectrons/150.

            colorIdx = 0
            for trackID,waveform in trackDict.items():
                label   = str('SimChannel ide Channel: %d, id: %d' % (simChannel.Channel(),trackID))
                print("--> label:",label)
                fig.add_trace(go.Scatter(x=np.arange(0,4096),y=waveform,name=label,line_shape="hvh",line=dict(color=colors[colorIdx],width=1)),row=wireID.Plane+1,col=1)
                colorIdx = (colorIdx + 1) % len(colors)
            
            break;
            
            #tdcIdeMap = simChannel.TDCIDEMap()
            #
            #for tdcIde in tdcIdeMap:
            #    #print("For TDCIDE tdc val:",tdcIde.first)
            #    numElectrons = 0
            #    trackID = -1
            #    for ide in tdcIde.second:
            #        numElectrons += ide.numElectrons
            #        if trackID < 0:
            #            trackID = ide.origTrackID
            #            
            #        #print("   --> numElectrons:",ide.numElectrons,", track id:",ide.trackID,", origTrackID:",ide.origTrackID)
            #        
            #if trackID in trackIDToMCTrack:
            #    start4Mom = trackIDToMCTrack[trackID].Start().Momentum()
            #    trackCosX = start4Mom.Px() / np.sqrt(start4Mom.Px()*start4Mom.Px()+start4Mom.Py()*start4Mom.Py()+start4Mom.Pz()*start4Mom.Pz())

    return
