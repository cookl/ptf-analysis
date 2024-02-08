#!/usr/bin/env python
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
# matplotlib.use('TkAgg')

print("Import matplotlib")
import scipy.stats as stats
import uproot3 as uproot
import numpy as np
import math 



def getTH1F(fileName, graphName):
    file = uproot.open(fileName)
    freq = np.array(file[graphName].values)
    bins = np.array(file[graphName].edges)
    bins = (bins[1:]+bins[:-1])/2.0
    return bins, freq

def getTH2F(fileName, graphName):
    file = uproot.open(fileName)
    values = np.array(file[graphName].values)
    # print(file[graphName].edges)
    xbins, ybins = np.array(file[graphName].edges,dtype = object)
    xbins = (xbins[1:]+xbins[:-1])/2.0
    ybins = (ybins[1:]+ybins[:-1])/2.0

    return xbins, ybins, values

def getPMTCoords(xbinCentres, ybinCentres, values):
    valueThreshold = 8
    values[values<valueThreshold] = 0
    
    #sum along each axis
    values_x = values.sum(axis=1)
    values_y = values.sum(axis=0)
    
    meanX = np.average(xbinCentres, weights =values_x)
    meanY = np.average(ybinCentres, weights =values_y)
    
    return meanX,meanY

def getXYCoordPMTPosition(pos):
    pmtPositionArray = np.array([[0.44746439,0.20756987],[0.45201267,0.30725607],[0.36634413,0.25269358],[0.36318949,0.16297828],[0.44738879,0.11218644],[0.53094878,0.15698004],[0.53155534,0.25496038],[0.45061428,0.3963719],[0.35594699,0.36829472],[0.29006354,0.30136363],[0.26234593,0.20776005],[0.29398638,0.11559652],[0.35488732,0.04815216],[0.45238402,0.02182653],[0.5428919,0.04878102],[0.60687074,0.11148197],[0.63501597,0.21287161],[0.61255901,0.30104175],[0.54385138,0.36338586]])
    return pmtPositionArray[pos]

def makePMTPosCircle(pos):
    x_centre, y_centre = getXYCoordPMTPosition(pos)
    bins = 1000
    x = np.zeros(bins,dtype=float)
    y = np.zeros(bins,dtype=float)
    
    r = 0.04
    for i in range(bins):
        theta = i*math.pi*2.0/bins
        x[i] = x_centre+(r*math.cos(theta))
        y[i] = y_centre+(r*math.sin(theta))
        # print("theta",theta,"x",x[i])
    
    # print("x,y")
    return x, y

def makeFeedInCircle(x_centre, y_centre):
    bins = 1000
    x = np.zeros(bins,dtype=float)
    y = np.zeros(bins,dtype=float)
    
    r = 0.02
    for i in range(bins):
        theta = i*math.pi*2.0/bins
        x[i] = x_centre+(r*math.cos(theta))
        y[i] = y_centre+(r*math.sin(theta))
        # print("theta",theta,"x",x[i])
    
    # print("x,y")
    return x, y
        
        

def getPositionFromChannel(ch, xPMT, yPMT):
    #this takes the value of x and y for the centre of each PMT for each individual position pmtPositionArray[position_i] = x,y
    #calibrated by run 
    
    minDistance =1E5
    min_pos = -1
    for i in range(19):
        pmtPos = getXYCoordPMTPosition(i)
        d = (pmtPos[0]-xPMT)**2 + (pmtPos[1]-yPMT)**2
        if(d<minDistance):
            min_pos=i
            minDistance =d
    
    if(minDistance>0.05):
        print("Error Channel ",ch, "is ",minDistance,"from closest PMT location")
        raise Exception("Error Channel ",ch, "is ",minDistance,"from closest PMT location")
    return min_pos

scanFileName = sys.argv[1]

#make the 2d scan plot
graphName = "SummedEfficiencyPlot"
xbins, ybins, values = getTH2F(scanFileName,graphName)
plt.pcolor(xbins, ybins, values.T)
plt.title(scanFileName + " Summed Efficiency")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.colorbar()
plt.savefig("TotalEff.png")
plt.close()

channelPMTCoords = np.zeros((19,2), dtype = float)
positionChannel = np.zeros((19), dtype = int)
# positionList = [17,6,0,18,1,7,8,2,9,3,10,11,4,12,13,5,14,15,16]

for i in range(19):
    #iterate over each channel
    graphName = "EfficiencyPlotCh"+str(i)
    xbins, ybins, values = getTH2F(scanFileName,graphName)
   
    #get the coords of the PMT from the channel
    xPMT, yPMT = getPMTCoords(xbins, ybins, values)
    channelPMTCoords[i,0]=xPMT
    channelPMTCoords[i,1]=yPMT
    #associate those coordinates with a PMT position
    pos = getPositionFromChannel(i, xPMT,yPMT)
    positionChannel[pos] = i
    
    xbins, ybins, values = getTH2F(scanFileName,graphName)

    
    plt.pcolor(xbins, ybins, values.T)
    plt.plot(xPMT,yPMT, marker ="X", color ="red")
    plt.title("Channel "+str(i))
    plt.savefig("Map_Channel_"+str(i)+".png")
    # plt.show(block=False)
    # input()
    plt.close()


fig = plt.subplots
for i in range(19):
    #iterate over each position
    x,y = makePMTPosCircle(i)
    # print(x)
    plt.plot(x,y, color = "black")
    
    #channel at position
    channel = positionChannel[i]
    xPMT,yPMT = channelPMTCoords[channel]
    print("xPMT,yPMT",xPMT,yPMT)
    plt.plot(xPMT,yPMT, marker ="X", color ="red")
    plt.annotate("Ch "+str(channel),(xPMT,yPMT),xytext=(-12, -12), textcoords='offset points')
    
    print("Position ",i, "is Channel", positionChannel[i])

x,y = makeFeedInCircle(0.57,0.2)
plt.plot(x,y, color = "grey", alpha = 0.5, label = "Feed through")

plt.legend()    
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.savefig("pmtMapping.png")
plt.close()


############## now do the darkRate + LED histograms histograms

if(len(sys.argv)>2):
    name = ["DarkRate","LED1","LED2","LED3"]
    for file in range(4):
        fileName = sys.argv[2+file]
        mean_pulseHeight = np.zeros(19, dtype =float)
        for i in range(19):
            #iterate over each channel
            graphName = "pulse_height_ch"+str(i)
            x, values = getTH1F(fileName, graphName)
            
            mean = np.average(x[x>15], weights = values[x>15])
            std = np.average((x[x>15]-mean)**2, weights = values[x>15] )
            std = math.sqrt(std)
            plt.step(x,values, where = "mid", label = "Data")
            mean_pulseHeight[i]=mean
            x_min = 0
            x_max = x[len(x)-1]
            x = np.linspace(x_min,x_max,num =1000)
            y = stats.norm.pdf(x, mean, std)*np.sum(values)
            # print("x",x)
            # print("y",y)
            plt.plot(x,y,label = "Normal fit")
            
            cumulative_sum = np.array([values[0:j].sum() for j in range(len(values))])/values.sum()
            x_lim_max = x[np.argmax(cumulative_sum>0.9)]
            plt.xlim((0,x_lim_max))
            plt.legend()
            plt.title(name[file] + " Channel " + str(i))
            plt.xlabel("Pulse Height/ mV")
            plt.ylabel("Freq.")
            plt.savefig(name[file]+"_"+str(i)+".png")
            plt.close()
            

        plt.scatter(channelPMTCoords[:,0],channelPMTCoords[:,1],c = mean_pulseHeight, s=100)
        plt.colorbar()
        x,y = makeFeedInCircle(0.57,0.2)
        plt.plot(x,y, color = "grey", alpha = 0.5, label = "Feed through")

        plt.title(name[file]+"Mean Pulse Height")
        plt.legend()
        plt.xlabel("x [m]")
        plt.ylabel("y [m]")
        plt.savefig(name[file]+"X_YScatter.png")
        plt.close()

    

