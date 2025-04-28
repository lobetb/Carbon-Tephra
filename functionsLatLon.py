# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 18:42:15 2020
Last modified on Mar 25 2024

@author: Benjamin Lobet
"""
#from joblib import Parallel, delayed
import utm
import numpy as np
import pandas as pd
import random as rd
import time
import scipy.io
from collections import Counter
import math
from os import path, mkdir,listdir,remove
from pathlib import Path
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
       
    data = data.loc[(data['VEI'] > 3)]
    data = data.loc[(data['VEI'] < 7)]
    
    data = data.reset_index(drop=True)

    return data

def get_ref_zone(data,refVolcano):
    indexLoc = data[data['Volcano Name'] == refVolcano].index.tolist()[0]
    coords = utm.from_latlon(data.loc[indexLoc,'Latitude'], data.loc[indexLoc,'Longitude'])
    numZone = coords[2]
    letterZone = coords[3]
    refLon = (math.trunc(coords[0]/1000))*1000
    refLat = (math.trunc(coords[1]/1000))*1000
    return (numZone,letterZone, refLat, refLon)


def get_ref_LatLon(data,refVolcano):
    indexLoc = data[data['Volcano Name'] == refVolcano].index.tolist()[0]
    refLat = data.loc[indexLoc, 'Latitude']
    refLon = data.loc[indexLoc, 'Longitude']
    return(refLat,refLon)


def get_prob(data, startYear, stopYear,thresholdYear,timeStepLength, mode): #Compute the probabilities of eruption (VEI 3 & 4) for each volcano
    print("Calculating probabilities for eruptions for each volcano...")
    counter0 = time.perf_counter()
    dataVEI4 = data.loc[(data['VEI'] == 4) & (data['Start Year'] >= thresholdYear)]
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
    counter1 = time.perf_counter()
    print("Done : " + str(round(counter1-counter0, 2)) + " seconds \n")
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

def get_stoch_eruptions(data, probabilities, startYear, stopYear,threshYear,refVolcano, mode):
    
    print("Computing eruptions according to probabilities and/or historical data...")
    counter0 = time.perf_counter()
    stopMonth = stopYear * 12
    startMonth = startYear * 12
    eruptions = pd.DataFrame(columns=("Volcano","Year","VEI"))    
    k=0
    
    if mode == "mixed":
        
        for j in range(len(probabilities)):
            x = rd.choices([1,0], [probabilities.loc[1,j],1-probabilities.loc[1,j]], k=(stopMonth-startMonth))
            y = get_index_positions(x,1)
            z = [probabilities.loc[3,j]]*len(y)
            x = [probabilities.loc[0,j]]*len(z)
    
            
            
            for i in range(len(x)):
                eruptions.loc[k] = [x[i],math.floor(y[i]/12)+startYear,z[i]]
                k+=1
                
        dataVEI56 = data.loc[(data['VEI'] > 4) & (data['VEI'] < 7)]
        for i in range(len(dataVEI56)):
    
        
            eruptions.loc[k] = [dataVEI56.iloc[i,1], dataVEI56.iloc[i,8], 
                                int(dataVEI56.iloc[i,5])]
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
    
            
                        
                        eruptions.loc[k] = [temp.iloc[l,1], temp.iloc[l,8], 
                                int(temp.iloc[l,5])]
                        k+=1
                        
    elif mode == "stochastic":
        
        for j in range(len(probabilities)):
            x = rd.choices([1,0], [probabilities.loc[1,j],1-probabilities.loc[1,j]], k=(stopMonth-startMonth))
            y = get_index_positions(x,1)
            z = [probabilities.loc[3,j]]*len(y)
            x = [probabilities.loc[0,j]]*len(z)
           
    
            
            for i in range(len(x)):
                eruptions.loc[k] = [x[i],math.floor(y[i]/12)+startYear,z[i]]
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
               probaEruptVEI4 += probabilities[i][1]
               listVEI4.append(probabilities[i])
           if probabilities[i][3] == 5:
               probaEruptVEI5 += probabilities[i][1]
               listVEI5.append(probabilities[i])
           if probabilities[i][3] == 6:
               probaEruptVEI6 += probabilities[i][1]
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

            
            eruptions.loc[k] = [volcano[0], math.floor(y[i]/12)+startYear,VEI[0]]
            k+=1
        
    eruptions = eruptions.sort_values(by=['Year'], ascending=True)
    eruptions = eruptions.reset_index(drop=True)
    
    counter1 = time.perf_counter()
    print("Done : " + str(round(counter1-counter0, 2)) + " seconds \n")
    return eruptions

    
def create_vei_files(outputFolder, refVolcano, data, refVEILatLonRel):
    print("Generating a VEI coordinates file for each eruption/volcano couple in UTM coordinates, if the file doesn't exist...")
    counter0 = time.perf_counter()
    savePath = outputFolder + "VEIs"
    VEI4 = refVEILatLonRel[0].copy()
    VEI5 = refVEILatLonRel[1].copy()
    VEI6 = refVEILatLonRel[2].copy()
    
    if not path.exists(savePath):
        mkdir(savePath)
    for i in range(len(data)):
        if not path.exists(savePath / Path(data.loc[i,'Volcano Name'] + "VEI" + str(int(data.loc[i,'VEI'])) + ".csv")):

            volcLat = data.loc[i,'Latitude']
            volcLon = data.loc[i,'Longitude']
            
            if(data.loc[i,'VEI']) == 4:
                dfCopy = VEI4.copy()
            elif(data.loc[i,'VEI']) == 5:
                dfCopy = VEI5.copy()
            elif(data.loc[i,'VEI']) == 6:
                dfCopy = VEI6.copy()

            dfUTM = pd.DataFrame(columns=("Easting","Northing","zNum","zLet","Probability"))
            
            for j in range(len(dfCopy)):
                try:
                    dfCopy.loc[j,'Latitude'] += volcLat
                    dfCopy.loc[j,'Longitude'] += volcLon
                    tempUTM = utm.from_latlon(dfCopy.loc[j,'Latitude'], dfCopy.loc[j,'Longitude'])
                    easting = (math.trunc(tempUTM[0]/1000))*1000
                    northing = (math.trunc(tempUTM[1]/1000))*1000
                    zNum = tempUTM[2]
                    zLet = tempUTM[3]
                    Prob = dfCopy.loc[j,'Probability']
                    
                    dfUTM.loc[len(dfUTM)] = [easting, northing, zNum, zLet, Prob]
                except: 
                    print(i, j, dfCopy.loc[j,'Longitude'], volcLon)
                    return
            del dfCopy
            dfUTM.to_csv(savePath + "/" + data.loc[i,'Volcano Name'] + 'VEI' + str(int(data.loc[i,'VEI'])) + '.csv', sep=',')
            print("Generated file for VEI " + str(int(data.loc[i,'VEI'])) + " for the volcano " + data.loc[i,'Volcano Name'])

    fileList = []
    uniqueID = (data["Volcano Name"] + ":" + data["VEI"].astype(int).astype(str)).unique()
    for i in range(len(uniqueID)):
        temp = uniqueID[i].split(":")
        fileList.append(temp[0] + "VEI" + temp[1] + ".csv")
    counter1 = time.perf_counter()
    print("Done : " + str(round(counter1-counter0, 2)) + " seconds \n")
    return fileList

def chk_file(fileName):
    if fileName=="invalid":
        print("You have made an invalid choice, try again and input only Y or N \n")
    a=input("A file with " + fileName + " has been found, do you want to load it ? This will take less time, but it is recommended if your list of volcanoes has changed. Load it ? (Y/N) \n > ")
    if a == "Y":
        return "Y"
    if a == "N":
        return "N"
    else:
        return chk_file("invalid")

def get_erupt_by_volc(eruptions, outputFolder):
    
    eruptListByVolc = []
    listNames = []
    for i in range(len(eruptions)):
        concName = eruptions.loc[i,'Volcano'] + str(eruptions.loc[i,'VEI'])
        listNames.append(concName)

    for i in np.unique(listNames):
        tempList = [str(i)]
        for j in range(len(listNames)):
            if listNames[j] == i:
                tempList.append(eruptions.loc[j,'Year'])
        eruptListByVolc.append(tempList)
    return eruptListByVolc

def get_result(result):
    global bigList
    bigList.append(result)
    
    
def add_eruptions_years(i,savePath,fileList, eruptions, startYear, eruptListByVolc):
    
    tempListAgr = []
    fileName = savePath + fileList[i]
    df = pd.read_csv(fileName, index_col=0)
    concName = str(fileList[i].rsplit('.',1)[0].split('VEI')[0]) + str(fileList[i].rsplit('.',1)[0].split('VEI')[1])
    for j in range(len(eruptListByVolc)):
        
        if eruptListByVolc.loc[0,j] == concName:
            for k in range(len(df)):
                tempList = []
                tempList.append([df.loc[k,'Northing'], df.loc[k,'Easting'], df.loc[k,'zNum'], df.loc[k,'zLet']])
                tempList[len(tempList)-1].append(eruptListByVolc.loc[1:len(eruptListByVolc[j]),j])
                tempListAgr.append(tempList)
    
    return tempListAgr

def add_eruptions_years_serial(outputFolder,fileList, eruptions, startYear, eruptListByVolc):
    print("Adding the eruptions years.")
    counter0 = time.perf_counter()
    bigList = []
    savePath = outputFolder + "VEIs/"
    for i in range(len(fileList)):
        fileName = savePath + fileList[i]
        df = pd.read_csv(fileName, index_col=0)
        concName = str(fileList[i].rsplit('.',1)[0].split('VEI')[0]) + str(fileList[i].rsplit('.',1)[0].split('VEI')[1])
        for j in range(len(eruptListByVolc)):
            
            if eruptListByVolc[j][0] == concName:
                for k in range(len(df)):
                    tempList = []
                    tempList.append([df.loc[k,'Northing'], df.loc[k,'Easting'], df.loc[k,'zNum'], df.loc[k,'zLet']])
                    tempList[len(tempList)-1].append(eruptListByVolc[j][1:len(eruptListByVolc[j])])
                    bigList.append(tempList)
    
    
    
    seen=set()
    bigListUnique = []
    for i in range(len(bigList)):
        conc = str(bigList[i][0][0]) + str(bigList[i][0][1]) + str(bigList[i][0][2]) + str(bigList[i][0][3])
        if conc not in seen:
            bigListUnique.append(bigList[i])
            bigListUnique[len(bigListUnique)-1][0][4].append(startYear)
            seen.add(conc)
        else:
            for j in range(len(bigListUnique)):
                if bigListUnique[j][0][0] == bigList[i][0][0] and bigListUnique[j][0][1] == bigList[i][0][1] and bigListUnique[j][0][2] == bigList[i][0][2] and bigListUnique[j][0][3] == bigList[i][0][3]:
                    for k in range(len(bigList[j][0][4])):
                        bigListUnique[j][0][4].append(bigList[j][0][4][k])
                bigListUnique[j][0][4].sort()                  

    counter1 = time.perf_counter()
    print("Done : " + str(round(counter1-counter0, 2)) + " seconds \n")
    return bigListUnique

def truncate_float(float_number, decimal_places):
    multiplier = 10 ** decimal_places
    return int(float_number * multiplier) / multiplier

def carbonAccumulation(x):
    res = 501.4*(x**(-0.55))
    return res

def compute_carbon(eruptionYears, startYear, stopYear, surfaceC, outputFolder, cellSize, carbonReduction):
    print("Computing the amount of carbon trapped by tephra...")
    counter0 = time.perf_counter()
    listC = []
    surfaceCList = [0]*len(eruptionYears)
    
    logC = [0]*(stopYear - startYear)
    for i in range(len(eruptionYears)):
        sumC = 0
        for j in range(1,len(eruptionYears[i][0][4])):
            timeDif = eruptionYears[i][0][4][j] - eruptionYears[i][0][4][j-1]
            if timeDif > 0:
                amountC, error = quad(carbonAccumulation, 0, timeDif)
                amountC = amountC * cellSize**2
                sumC += amountC*carbonReduction
                logC[eruptionYears[i][0][4][j]] += sumC
        listC.append(sumC)
        
    if surfaceC == "yes":
        
        for i in range(len(eruptionYears)):
            timeDif = stopYear - eruptionYears[i][0][4][-1]
            if timeDif > 0:
                amountC, error = quad(carbonAccumulation, 0, timeDif)
                amountC = amountC * cellSize**2
                surfaceCList[i] = amountC
                
    counter1 = time.perf_counter()
    print("Done : " + str(round(counter1-counter0, 2)) + " seconds \n")
    return listC, logC, surfaceCList

def create_dataset_latlon(eruptionYears, carbonList, surfaceCList):
    print("Compiling the dataset...")
    counter0 = time.perf_counter()
    df = pd.DataFrame(columns=('Latitude','Longitude','Trapped C','Surface C'))
    
    for i in range(len(eruptionYears)):
        (lat,lon) = utm.to_latlon(eruptionYears[i][0][1], eruptionYears[i][0][0], eruptionYears[i][0][2],eruptionYears[i][0][3])
        df.loc[i] = [lat, lon, carbonList[i], surfaceCList[i]]
    counter1 = time.perf_counter()
    print("Done : " + str(round(counter1-counter0, 2)) + " seconds \n")   
    return df
    
def save_results(outputFolder, dataSetLatLon, logC, eruptions):
    print("Saving the results to file...")
    counter0 = time.perf_counter()
    if not path.exists(outputFolder + "Runs/") :
        mkdir(outputFolder + "Runs/")
    dataSetLatLon.to_csv(outputFolder + "Runs/run" + str(len(listdir(outputFolder + "Runs/"))+1) + ".csv", index=False)
    if not path.exists(outputFolder + "LogC/"):
        mkdir(outputFolder + "LogC/")
    np.savetxt(outputFolder + "LogC/logCrun" + str(len(listdir(outputFolder + "LogC/"))+1) + ".csv",logC, delimiter=",")
    if not path.exists(outputFolder + "LogErupt/"):
        mkdir(outputFolder + "LogErupt/")
    eruptions.to_csv(outputFolder + "LogErupt/eruptRun" + str(len(listdir(outputFolder + "LogErupt/"))+1) + ".csv", index=False)  
    count = len(listdir(outputFolder + "LogC/"))+1
    counter1 = time.perf_counter()
    print("Done : " + str(round(counter1-counter0, 2)) + " seconds \n") 
    
    return count

def import_data(inputFileFolder,inputFile, VEI4MAT,limits,refVolcano):
    counter0 = time.perf_counter()
    print("Importing data... \n")
    fullPath = inputFileFolder + inputFile

    data = pd.read_excel(fullPath, sheet_name="Eruption List")

    mat = scipy.io.loadmat(inputFileFolder + VEI4MAT)

    atacazoVEI4 = np.array(mat['atacazo_vei4'])
    atacazoVEI5 = np.array(mat['atacazo_vei5'])
    atacazoVEI6 = np.array(mat['atacazo_vei6'])
    refVEI = [atacazoVEI4,atacazoVEI5,atacazoVEI6]
        
   
    data = apply_coord_constraints(data, limits)
    numZone,letterZone, refLat, refLon = get_ref_zone(data, refVolcano)
    counter1 = time.perf_counter()
    print("Done : " + str(round(counter1-counter0, 2)) + " seconds \n")
    return refVEI, data, refLat, refLon, numZone, letterZone
    
def convertLatLon(data,refVolcano,refVEI, probThreshold):
    counter0 = time.perf_counter()
    print("Converting tephra data in latlon coordinates...")
    indexLoc = data[data['Volcano Name'] == refVolcano].index.tolist()[0]
    refLat = data.loc[indexLoc,'Latitude'] 
    refLon = data.loc[indexLoc,'Longitude']
    atacazoVEI4 = refVEI[0]
    atacazoVEI5 = refVEI[1]
    atacazoVEI6 = refVEI[2]
    
    refZone = get_ref_zone(data, refVolcano)
    refZoneNumber = refZone[0]
    refZoneLetter = refZone[1]
    
    atacazoVEI4LatLon = pd.DataFrame(columns=("Latitude","Longitude","Probability"))
    k = 0
    for i in range(len(atacazoVEI4)):
        number = refZoneNumber
        letter = refZoneLetter
        northing = atacazoVEI4[i,1]
        easting = atacazoVEI4[i,0]

        if atacazoVEI4[i,2] > probThreshold:
            LatLon = utm.to_latlon(easting,northing,number,letter)
            atacazoVEI4LatLon.loc[k] = (LatLon[0], LatLon[1], atacazoVEI4[i,2])
            k += 1
            
    for i in range(len(atacazoVEI4LatLon)):
        atacazoVEI4LatLon.loc[i,'Latitude'] -= refLat
        atacazoVEI4LatLon.loc[i,'Longitude'] -= refLon
            
    atacazoVEI5LatLon = pd.DataFrame(columns=("Latitude","Longitude","Probability"))
    k = 0
    for i in range(len(atacazoVEI5)):
        number = refZoneNumber
        letter = refZoneLetter
        northing = atacazoVEI5[i,1]
        easting = atacazoVEI5[i,0]
        if atacazoVEI5[i,1] > 10000000:
            letter = chr(ord(letter) + 1)
            northing -= 10000000
        if atacazoVEI5[i,2] > probThreshold:
            LatLon = utm.to_latlon(easting,northing,number,letter)
            atacazoVEI5LatLon.loc[k] = (LatLon[0], LatLon[1], atacazoVEI5[i,2])
            k += 1
    for i in range(len(atacazoVEI5LatLon)):
        atacazoVEI5LatLon.loc[i,'Latitude'] -= refLat
        atacazoVEI5LatLon.loc[i,'Longitude'] -= refLon
    
    atacazoVEI6LatLon = pd.DataFrame(columns=("Latitude","Longitude","Probability"))
    k = 0
    for i in range(len(atacazoVEI6)):
        number = refZoneNumber
        letter = refZoneLetter
        northing = atacazoVEI6[i,1]
        easting = atacazoVEI6[i,0]
        if atacazoVEI6[i,1] > 10000000:
            letter = chr(ord(letter) + 1)
            northing -= 10000000
        if atacazoVEI6[i,2] > probThreshold:
            LatLon = utm.to_latlon(easting,northing,number,letter)
            atacazoVEI6LatLon.loc[k] = (LatLon[0], LatLon[1], atacazoVEI6[i,2])
            k += 1
    for i in range(len(atacazoVEI6LatLon)):
        atacazoVEI6LatLon.loc[i,'Latitude'] -= refLat
        atacazoVEI6LatLon.loc[i,'Longitude'] -= refLon

    refVEILatLonRel = [atacazoVEI4LatLon,atacazoVEI5LatLon,atacazoVEI6LatLon]
    counter1 = time.perf_counter()
    print("Done : " + str(round(counter1-counter0, 2)) + " seconds \n")
    return refVEILatLonRel

def write_parameters(outputFolder, timeStepLength, inputFile, VEI4MAT, refVolcano, carbonReduction, probThreshold, 
                   cellSize, thresholdYear, stopYear, startYear, limits, surfaceC, mode):
    counter0 = time.perf_counter()
    print("Writing parameters to file...")
    if not path.exists(outputFolder):
        mkdir(outputFolder)
    if path.exists(outputFolder + "parameters.txt"):
        remove(outputFolder + "parameters.txt")
    paramLog = open(outputFolder + "parameters.txt", 'a')
    paramLog.write("timeStepLength = " + str(timeStepLength) + "\n")
    paramLog.write("inputFile = " + str(inputFile) + "\n")
    paramLog.write("VEI4MAT = " + str(VEI4MAT) + "\n")
    paramLog.write("refVolcano = " + str(refVolcano) + "\n")
    paramLog.write("carbonReduction = " + str(carbonReduction) + "\n")
    paramLog.write("probThreshold = " + str(probThreshold) + "\n")
    paramLog.write("cellSize = " + str(cellSize) + "\n")
    paramLog.write("thresholdYear = " + str(thresholdYear) + "\n")
    paramLog.write("stopYear = " + str(stopYear) + "\n")
    paramLog.write("startYear = " + str(startYear) + "\n")
    paramLog.write("limits = " + str(limits) + "\n")
    paramLog.write("surfaceC = " + str(surfaceC) + "\n")
    paramLog.write("mode = " + str(mode) + "\n")
    paramLog.close()
    counter1 = time.perf_counter()
    print("Done : " + str(round(counter1-counter0, 2)) + " seconds \n")