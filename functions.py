# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 18:42:15 2020

@author: Ben
"""
import utm
import numpy as np
import pandas as pd
import random as rd
from collections import Counter
import math
from os import path, mkdir,listdir
from pathlib import Path
from os.path import isfile, join
from scipy.integrate import quad


def get_index_positions(list_of_elems, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
    index_pos_list = []
    index_pos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            index_pos = list_of_elems.index(element, index_pos)
            # Add the index position in list
            index_pos_list.append(index_pos)
            index_pos += 1
        except ValueError as e:
            break
    return index_pos_list

def apply_coord_constraints(data, coordConstraints):
    if not coordConstraints[0] == None:
        data = data.loc[(data['Latitude'] < coordConstraints[0])]
    if not coordConstraints[1] == None:
        data = data.loc[(data['Latitude'] > coordConstraints[1])]
    if not coordConstraints[2] == None:
        data = data.loc[(data['Longitude'] < coordConstraints[2])]
    if not coordConstraints[3] == None:
        data = data.loc[(data['Longitude'] > coordConstraints[3])]
    data = data.reset_index(drop=True)

    return data

def get_ref_zone(data,refVolcano):
    indexLoc = data[data['Volcano Name'] == refVolcano].index.tolist()[0]
    coords = utm.from_latlon(data.iloc[indexLoc,22], data.iloc[indexLoc,23])
    numZone = coords[2]
    letterZone = coords[3]
    refLon = coords[0]
    refLat = coords[1]
    return (numZone,letterZone, refLat, refLon)


def get_prob(data, startYear, stopYear,thresholdYear,timeStepLength, mode): #Compute the probabilities of eruption (VEI 3 & 4) for each volcano
  
    dataVEI4 = data.loc[(data['VEI'] > 3) & (data['VEI'] < 5) & (data['Start Year'] >= thresholdYear)]
    counterVolcanoesVEI4 = Counter(dataVEI4['Volcano Name'])
    probEruptions = []
    
    timeStepsProbs = (stopYear - thresholdYear) * (12 * timeStepLength)
    
    for i in counterVolcanoesVEI4.keys():
        newLine = [i, counterVolcanoesVEI4[i]/timeStepsProbs,
                                len(dataVEI4.loc[(dataVEI4['VEI'] == 4) 
                                                  & (dataVEI4['Volcano Name'] == i)]), 4]
        probEruptions.append(newLine)
    
    if not mode == "mixed":
        
        timeStepsProbsVEI56 = (stopYear - startYear) * (12 *timeStepLength)
        dataVEI5 = data.loc[(data['VEI'] == 5)]
        counterVolcanoesVEI5 = Counter(dataVEI5['Volcano Name'])
        
        
        for i in counterVolcanoesVEI5.keys():
            newLine = [i, counterVolcanoesVEI5[i]/timeStepsProbsVEI56, len(dataVEI5.loc[(dataVEI5['VEI'] == 5) & (dataVEI5['Volcano Name'] == i)]), 5]
            probEruptions.append(newLine)
            
        dataVEI6 = data.loc[(data['VEI'] == 6)]
        counterVolcanoesVEI6 = Counter(dataVEI6['Volcano Name'])
        
        
        for i in counterVolcanoesVEI6.keys():
            newLine = [i, counterVolcanoesVEI6[i]/timeStepsProbsVEI56, len(dataVEI6.loc[(dataVEI6['VEI'] == 6) & (dataVEI6['Volcano Name'] == i)]), 6]
            probEruptions.append(newLine)
            
        dataVEI4 = data.loc[(data['VEI'] == 4)]
        counterVolcanoesVEI4 = Counter(dataVEI4['Volcano Name'])
        for i in counterVolcanoesVEI4.keys():
            if len(dataVEI4.loc[(dataVEI4['Start Year'] >= thresholdYear) & 
                                (dataVEI4['Volcano Name'] == i)]) == 0:
                if len(dataVEI4.loc[(dataVEI4['Start Year'] < thresholdYear) & 
                                (dataVEI4['Volcano Name'] == i)]) > 0:
                    temp = dataVEI4.loc[(dataVEI4['Start Year'] < thresholdYear) & 
                                (dataVEI4['Volcano Name'] == i)]
                    newLine = [i, counterVolcanoesVEI4[i]/timeStepsProbsVEI56, len(temp), 4]
                    probEruptions.append(newLine)
        
    return probEruptions

def get_utm_zone(lat,long):
    x1 = long + 180
    x2 = x1 / 6
    numZone = round(x2)
    NS = None
    if lat >= 0:
        NS = "N"
    elif lat < 0:
        NS = "S"
    
    return numZone, NS

def get_stoch_eruptions(data, probabilities, startYear, stopYear,threshYear,refZone, mode):
    
    stopMonth = stopYear * 12
    startMonth = startYear * 12
    eruptions = pd.DataFrame(columns=("Volcano","Year","VEI", "Lat","Lon"))    
    k=0
    
    if mode == "mixed":
        
        for j in range(len(probabilities)):
            x = rd.choices([1,0], [probabilities[j][1],1-probabilities[j][1]], k=(stopMonth-startMonth))
            y = get_index_positions(x,1)
            z = [probabilities[j][3]]*len(y)
            x = [probabilities[j][0]]*len(z)
            locIndex = data.index[data['Volcano Name'] == x[0]].tolist()[0]
    
            coords = utm.from_latlon(data.iloc[locIndex,22],data.iloc[locIndex,23],
                                     force_zone_number = refZone[0],force_zone_letter = refZone[1])
            latC = coords[1]
            lonC = coords[0]
            if data.iloc[locIndex,22] >= 0:
                latC = latC+10000000
            lat = [latC]*len(z)
            lon = [lonC]*len(z)
            
            for i in range(len(x)):
                eruptions.loc[k] = [x[i],math.floor(y[i]/12)+startYear,z[i], lat[i], lon[i]]
                k+=1
                
        dataVEI56 = data.loc[(data['VEI'] > 4) & (data['VEI'] < 7)]
        for i in range(len(dataVEI56)):
    
            coords = utm.from_latlon(dataVEI56.iloc[i,22],dataVEI56.iloc[i,23],
                                     force_zone_number = refZone[0],force_zone_letter = refZone[1])
            latC = coords[1]
            lonC = coords[0]
            if dataVEI56.iloc[i,22] >= 0:
                latC = latC+10000000
            eruptions.loc[k] = [dataVEI56.iloc[i,1], dataVEI56.iloc[i,8], 
                                int(dataVEI56.iloc[i,5]), latC, lonC]
            k+=1
        dataVEI4 = data.loc[(data['VEI'] == 4)]
        counterVolcanoes = Counter(dataVEI4['Volcano Name'])
        for i in counterVolcanoes.keys():
            if len(dataVEI4.loc[(dataVEI4['Start Year'] >= threshYear) & 
                                (dataVEI4['Volcano Name'] == i)]) == 0:
                if len(dataVEI4.loc[(dataVEI4['Start Year'] < threshYear) & 
                                (dataVEI4['Volcano Name'] == i)]) > 0:
                    temp = dataVEI4.loc[(dataVEI4['Start Year'] < threshYear) & 
                                (dataVEI4['Volcano Name'] == i)]
                    for l in range(len(temp)):
    
                        coords = utm.from_latlon(temp.iloc[l,22],temp.iloc[l,23],
                                     force_zone_number = refZone[0],force_zone_letter = refZone[1])
                        latC = coords[1]
                        lonC = coords[0]
                        if temp.iloc[l,22] >= 0:
                            latC = latC+10000000
                        
                        eruptions.loc[k] = [temp.iloc[l,1], temp.iloc[l,8], 
                                int(temp.iloc[l,5]),latC,lonC]
                        k+=1
                        
    elif mode == "stochastic":
        
        for j in range(len(probabilities)):
            x = rd.choices([1,0], [probabilities[j][1],1-probabilities[j][1]], k=(stopMonth-startMonth))
            y = get_index_positions(x,1)
            z = [probabilities[j][3]]*len(y)
            x = [probabilities[j][0]]*len(z)
            locIndex = data.index[data['Volcano Name'] == probabilities[j][0]].tolist()[0]
    
            coords = utm.from_latlon(data.iloc[locIndex,22],data.iloc[locIndex,23],
                                     force_zone_number = refZone[0],force_zone_letter = refZone[1])
            latC = coords[1]
            lonC = coords[0]
            if data.iloc[locIndex,22] >= 0:
                latC = latC+10000000
            lat = [latC]*len(z)
            lon = [lonC]*len(z)
            
            for i in range(len(x)):
                eruptions.loc[k] = [x[i],math.floor(y[i]/12)+startYear,z[i], lat[i], lon[i]]
                k+=1
                       
                
    elif mode == "sequential":
        
        probaEruptTot = 0
        probaEruptVEI4 = 0
        probaEruptVEI5 = 0
        probaEruptVEI6 = 0
        listVEI4 = []
        listVEI5 = []
        listVEI6 = []
        
        for i in range(len(probabilities)):
            probaEruptTot += probabilities[i][1]
            if probabilities[i][3] == 4:
                probaEruptVEI4 += probabilities[i][3]
                listVEI4.append(probabilities[i])
            if probabilities[i][3] == 5:
                probaEruptVEI5 += probabilities[i][3]
                listVEI5.append(probabilities[i])
            if probabilities[i][3] == 6:
                probaEruptVEI6 += probabilities[i][3]
                listVEI6.append(probabilities[i])

        x = rd.choices([1,0],[probaEruptTot, 1-probaEruptTot], k=stopMonth-startMonth)
        y = get_index_positions(x,1)
        
        for i in range(len(y)):
            VEI = rd.choices([4,5,6], [probaEruptVEI4, probaEruptVEI5, probaEruptVEI6])

            vol = []
            prob = []
            if VEI == [4]:
                for j in range(len(listVEI4)):
                    vol.append(listVEI4[j][0])
                    prob.append(listVEI4[j][1])
            elif VEI == [5]:
                for j in range(len(listVEI5)):
                    vol.append(listVEI5[j][0])
                    prob.append(listVEI5[j][1])
            elif VEI == [6]:
                for j in range(len(listVEI6)):
                    vol.append(listVEI6[j][0])
                    prob.append(listVEI6[j][1])
                    
            volcano = rd.choices(vol, weights=prob)
            locIndex = data.index[data['Volcano Name'] == volcano[0]].tolist()[0]
            coords = utm.from_latlon(data.iloc[locIndex,22],data.iloc[locIndex,23],
                                     force_zone_number = refZone[0],force_zone_letter = refZone[1])
            latC = coords[1]
            lonC = coords[0]
            if data.iloc[locIndex,22] >= 0:
                latC = latC+10000000
            
            eruptions.loc[k] = [volcano[0], math.floor(y[i]/12)+startYear,VEI[0], latC, lonC]
            k+=1
        
    eruptions = eruptions.sort_values(by=['Year'], ascending=False)
    return(eruptions)

def create_vei_files(inputFileFolder, refVolcano, eruptions, refVEI, refZone):
    refLat = refZone[2]
    refLon = refZone[3]
    
    savePath = inputFileFolder + "VEIs"
    for i in range(len(eruptions)):
        
        if not path.exists(savePath):
            mkdir(savePath)
        if not path.exists(savePath / Path(eruptions.iloc[i,0] + "VEI" + str(eruptions.iloc[i,2]) + ".csv")):
            latDecal = eruptions.iloc[i,3] - refLat
            lonDecal = eruptions.iloc[i,4] - refLon
            if(eruptions.iloc[i,2]) == 4:
                matCopy = refVEI[0].copy()
            elif(eruptions.iloc[i,2]) == 5:
                matCopy = refVEI[1].copy()
            elif(eruptions.iloc[i,2]) == 6:
                matCopy = refVEI[2].copy()
            
            matCopy[:,1] += latDecal
            matCopy[:,0] += lonDecal
            
            np.savetxt((savePath / Path(eruptions.iloc[i,0] + "VEI" + str(eruptions.iloc[i,2]) + ".csv")),
                        matCopy, delimiter=",")
    fileList = []
    uniqueID = (eruptions["Volcano"] + "." + eruptions["VEI"].astype(str)).unique()
    for i in range(len(uniqueID)):
        temp = uniqueID[i].split(".")
        fileList.append(temp[0] + "VEI" + temp[1] + ".csv")
    return fileList
    
def create_grid(inputFileFolder, fileList, cellSize, outputFolder):
    if not path.exists(outputFolder):
        mkdir(outputFolder)
    if not path.exists(outputFolder + "grid.npy"):
        minLat = None
        maxLat = None
        minLon = None
        maxLon = None
        savePath = inputFileFolder + "VEIs"
        
        for i in range(len(fileList)):
            temp1 = pd.read_csv(savePath / Path(fileList[i]))
            temp = pd.DataFrame()
            temp = temp1[temp1.iloc[:,2] > 0.8]
            if minLat == None:
                minLat = min(temp.iloc[:,1])
            elif minLat > min(temp.iloc[:,1]):
                minLat = min(temp.iloc[:,1])
                
            if maxLat == None:
                maxLat = max(temp.iloc[:,1])
            elif maxLat < max(temp.iloc[:,1]):
                maxLat = max(temp.iloc[:,1])
            
            if minLon == None:
                minLon = min(temp.iloc[:,0])
            elif minLon > min(temp.iloc[:,0]):
                minLon = min(temp.iloc[:,0])
                
            if maxLon == None:
                maxLon = max(temp.iloc[:,0])
            elif maxLon < max(temp.iloc[:,0]):
                maxLon = max(temp.iloc[:,0])
                
        mapMatrix = np.empty((math.ceil(maxLat/cellSize) - math.floor(minLat/cellSize),
                                  math.ceil(maxLon/cellSize) - math.floor(minLon/cellSize)), dtype=object)
        for i in range(mapMatrix.shape[0]):
            for j in range(mapMatrix.shape[1]):
                mapMatrix[i,j] = []
                mapMatrix[i,j].append((math.ceil(minLat/cellSize)*cellSize)+(i*cellSize))
                mapMatrix[i,j].append((math.ceil(minLon/cellSize)*cellSize)+(j*cellSize))
        np.save(outputFolder + "grid", mapMatrix)
        coords = [minLat,minLon]
        np.save(outputFolder + "coords", coords)
        
    else:
        mapMatrix = np.load(outputFolder + "grid.npy", allow_pickle=True)
        coords = np.load(outputFolder + "coords.npy")
        minLat = coords[0]
        minLon = coords[1]
    
    return mapMatrix, math.ceil(minLat/cellSize)*cellSize, math.ceil(minLon/cellSize)*cellSize

def add_eruptions_to_grid(inputFileFolder,fileList, eruptions, grid, probThreshold, minLat, minLon, cellSize):
    veis = dict()
    savePath = inputFileFolder + "VEIs"
    for i in range(len(fileList)):
        veis[fileList[i].split("/")[-1].rsplit(".",1)[0]] = pd.read_csv(savePath / Path(fileList[i]))
        
    for i in range(len(eruptions)):
        vei = eruptions.iloc[i,0] + "VEI" + str(eruptions.iloc[i,2])
        values = veis[vei]
        valuesThresh = values[(values.iloc[:,2] >= probThreshold)]
        
        for j in range(len(valuesThresh)):
            grid[(math.floor((valuesThresh.iloc[j,1]-minLat)/cellSize)),math.floor((valuesThresh.iloc[j,0]-minLon)/cellSize)].append(eruptions.iloc[i,1])
        

def carbonAccumulation(x):
    res = 501.4*(x**(-0.55))
    return res

def get_carbon_grid(grid, startYear, stopYear):
    carbonGrid = np.empty((len(grid[:]),len(grid[0])))
    logC = [0] * (stopYear - startYear)
    for i in range(len(grid[:])):
        for j in range(len(grid[0])):
            sumC = 0
            if len(grid[i,j]) > 2:
                if not grid[i,j][-1] == startYear:
                    grid[i,j].append(startYear)
                k = len(grid[i,j]) -1
                while k > 2:
                    timeDif = grid[i,j][k-1] - grid[i,j][k]
                    if timeDif > 0:
                        amountC, error = quad(carbonAccumulation, 0, timeDif)
                        logC[grid[i,j][k-1]] += amountC/2
                        sumC += amountC/2
                    k -= 1
            carbonGrid[i,j] = sumC
            
    return carbonGrid, logC

def save_results(outputFolder, carbonGrid, logC, eruptions):
    if not path.exists(outputFolder + "Runs/") :
        mkdir(outputFolder + "Runs/")
    np.savetxt(outputFolder + "Runs/run" + str(len(listdir(outputFolder + "Runs/"))+1) + ".csv",carbonGrid, delimiter=",")
    if not path.exists(outputFolder + "LogC/"):
        mkdir(outputFolder + "LogC/")
    np.savetxt(outputFolder + "LogC/logCrun" + str(len(listdir(outputFolder + "LogC/"))+1) + ".csv",logC, delimiter=",")
    if not path.exists(outputFolder + "LogErupt/"):
        mkdir(outputFolder + "LogErupt/")
    eruptions[['Volcano','Year','VEI']].to_csv(outputFolder + "LogErupt/eruptRun" + str(len(listdir(outputFolder + "LogErupt/"))+1) + ".csv", header=False, index=False)  
    count = len(listdir(outputFolder + "LogC/"))+1
    
    return count
        