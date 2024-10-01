# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:34:06 2020

@author: Benjamin Lobet, benjaminlobet@gmail.com
"""

"""
Import necessary packages
"""
from os import chdir, path, mkdir, remove
chdir('g:/mon drive/Carbon-Tephra-master/')
import functionsLatLon as f2
import time

"""
Declare variables
"""
totalCountInit = time.perf_counter() #Initiate counter for whole simulation

timeStepLength = 1            #Time step for the simulation in months
inputFileFolder = "g:/mon drive/Carbon-Tephra-master/Input files/"
inputFile = "GVP_Eruption_Results.xls.xlsx"          #File containing infos for the eruptions, from GVP website
outputFolder = "g:/mon drive/Carbon-Tephra-master/Results/"
VEI4MAT = "atacazo_vei4.mat"
refVolcano = 'Atacazo'
carbonReduction = 0.5   #Proportion of carbon loss in buried soils
probThreshold = 0.8     #Threshold of probability for the isopach
cellSize = 1000         #Size of one side of each cell, in meters
thresholdYear = 1500    #Year above which the record is complete
stopYear = 2023        #Year when the records end
startYear = -10000      #Start year of the holocène
northLimit = None
southLimit = None
eastLimit = None
westLimit = None
limits = [northLimit, southLimit, eastLimit, westLimit] # Geographical limits (in lat/lon) to apply to the dataset to exclude some data.

surfaceC = "yes"        # Compute the accumulation of C in the surface soil. "yes" if it needs to be computed. "no" if it doesn't, or "fuck off",
                        # or anything except "yes" actually.
    
mode = "sequential"     # Mode of inference of the eruptions. "stochastic" does a full stochastic approach
                        # for VEI4-5-6, "mixed" uses the historical data for VEI 5 and 6 and stochastic for
                        # VEI4. 
                        # "sequential" takes a more sequential approach : first computing
                        # if there's an eruption, then the VEI, then attribution to a volcano. "sequential" uses the stochastic
                        # approach.

"""
Function calls
"""

#Import the data from the VEI reference file and the GVP table
refVEI, data, refLat, refLon, numZone, letterZone = f2.import_data(inputFileFolder,inputFile, VEI4MAT,limits,refVolcano)

#Convert coordinates in Lat/Lon coordinate system
refVEILatLonRel = f2.convertLatLon(data,refVolcano,refVEI, probThreshold)

#Write parameters of the simulation to file parameters.txt
f2.write_parameters(outputFolder, timeStepLength, inputFile, VEI4MAT, refVolcano, carbonReduction, probThreshold, 
                   cellSize, thresholdYear, stopYear, startYear, limits, surfaceC, mode)

#Compute the probabilities of eruption for each volcano/VEI combination
probabilities = f2.get_prob(data, startYear, stopYear, thresholdYear,timeStepLength, mode)

#Create a VEI file per volcano/VEI combination and return the list
fileList = f2.create_vei_files(outputFolder, refVolcano, data, refVEILatLonRel)

#Get the list of eruptions computed from the probabilities
eruptions = f2.get_stoch_eruptions(data, probabilities, startYear, stopYear, thresholdYear, refVolcano, mode)

#Create a list of eruptions per volcano
eruptListByVolc = f2.get_erupt_by_volc(eruptions, outputFolder)

#Creates a list containing, for each coordinate set, the years of eruptions
eruptionYears = f2.add_eruptions_years_serial(outputFolder,fileList, eruptions, startYear, eruptListByVolc)

#Computes, for each coordinates, the amount of C captured by the tephra + the C accumulated on the surface
carbonList, logC, surfaceCList = f2.compute_carbon(eruptionYears, startYear, stopYear, surfaceC, outputFolder, cellSize, carbonReduction)

#Create dataset for plotting a map (see below)
dataSetLatLon = f2.create_dataset_latlon(eruptionYears, carbonList, surfaceCList)

#save results
count = f2.save_results(outputFolder, dataSetLatLon, logC, eruptions)

totalCountFinal = time.perf_counter()

print("Whole computation in " + mode + " mode : " + str(round(totalCountFinal-totalCountInit, 2)) + " seconds \n")

        
"""
Plot on interactive map
"""
"""
import folium
import branca.colormap as cm
linear = cm.LinearColormap(["yellow", "orange", "red"], vmin=min(dataSetLatLon['Trapped C']), vmax=max(dataSetLatLon['Trapped C']))
linear2 = cm.LinearColormap(["yellow", "orange", "red"], vmin=min(dataSetLatLon['Surface C']), vmax=max(dataSetLatLon['Surface C']))
loc_center = [dataSetLatLon['Latitude'].mean(), dataSetLatLon['Longitude'].mean()]
map1 = folium.Map(location=loc_center, tiles='Openstreetmap', zoom_start=7, control_scale=True)

fg = folium.FeatureGroup(name="Trapped C", control=True).add_to(map1)
for index, loc in dataSetLatLon.iterrows():
    folium.Circle([loc['Latitude'],loc['Longitude']],radius=500,weight=5, stroke=False, fill_color=linear(loc['Trapped C']), 
                  popup=(loc['Trapped C']), fill_opacity=0.5).add_to(fg)
if surfaceC == "yes":
    fg2 = folium.FeatureGroup(name="Surface C", control=True).add_to(map1)
    for index, loc in dataSetLatLon.iterrows():
        folium.Circle([loc['Latitude'],loc['Longitude']],radius=500,weight=5, stroke=False, fill_color=linear2(loc['Surface C']), 
                      popup=(loc['Surface C']), fill_opacity=0.5).add_to(fg2)
folium.LayerControl().add_to(map1)
map1.save('map2.html')
"""


"""
Batch execution
""" 
"""
mode = "mixed"
for j in (0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9):
    probThreshold = j
    outputFolder = "h:/mon drive/Carbon-Tephra-master/Results/"+ str(mode) +"/" + str(j) + "/"
    
    
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






    for i in range(101):
        
        counter0 = time.perf_counter()
        refVEI, data, refLat, refLon, numZone, letterZone = f2.import_data(inputFileFolder,inputFile, VEI4MAT,limits,refVolcano)
        refVEILatLonRel = f2.convertLatLon(data,refVolcano,refVEI, probThreshold)
        f2.write_parameters(outputFolder, timeStepLength, inputFile, VEI4MAT, refVolcano, carbonReduction, probThreshold, 
                           cellSize, thresholdYear, stopYear, startYear, limits, surfaceC, mode)
        probabilities = f2.get_prob(data, startYear, stopYear, thresholdYear,timeStepLength, mode)
        fileList = f2.create_vei_files(outputFolder, refVolcano, data, refVEILatLonRel)
        eruptions = f2.get_stoch_eruptions(data, probabilities, startYear, stopYear, thresholdYear, refVolcano, mode)
        eruptListByVolc = f2.get_erupt_by_volc(eruptions, outputFolder)
        eruptionYears = f2.add_eruptions_years_serial(outputFolder,fileList, eruptions, startYear, eruptListByVolc)
        carbonList, logC, surfaceCList = f2.compute_carbon(eruptionYears, startYear, stopYear, surfaceC, outputFolder, cellSize, carbonReduction)

        
        counter1 = time.perf_counter()
        print(str(count) + ": " + str(counter1-counter0) + " secondes")

"""
"""
mode = "sequential"
for j in (0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9):
    probThreshold = j
    outputFolder = "h:/mon drive/Carbon-Tephra-master/Results/"+ str(mode) +"/" + str(j) + "/"
    
    
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






    for i in range(40):
        
        counter0 = time.perf_counter()
        refVEI, data, refLat, refLon, numZone, letterZone = f2.import_data(inputFileFolder,inputFile, VEI4MAT,limits,refVolcano)
        refVEILatLonRel = f2.convertLatLon(data,refVolcano,refVEI, probThreshold)
        f2.write_parameters(outputFolder, timeStepLength, inputFile, VEI4MAT, refVolcano, carbonReduction, probThreshold, 
                           cellSize, thresholdYear, stopYear, startYear, limits, surfaceC, mode)
        probabilities = f2.get_prob(data, startYear, stopYear, thresholdYear,timeStepLength, mode)
        fileList = f2.create_vei_files(outputFolder, refVolcano, data, refVEILatLonRel)
        eruptions = f2.get_stoch_eruptions(data, probabilities, startYear, stopYear, thresholdYear, refVolcano, mode)
        eruptListByVolc = f2.get_erupt_by_volc(eruptions, outputFolder)
        eruptionYears = f2.add_eruptions_years_serial(outputFolder,fileList, eruptions, startYear, eruptListByVolc)
        carbonList, logC, surfaceCList = f2.compute_carbon(eruptionYears, startYear, stopYear, surfaceC, outputFolder, cellSize, carbonReduction)

        
        counter1 = time.perf_counter()
        print(str(count) + ": " + str(counter1-counter0) + " secondes")
        """
