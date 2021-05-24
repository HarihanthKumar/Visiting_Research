import pandas as pd 
import numpy as np 
import os
import math
#load the two databases
cpgDataframe = pd.read_csv('/home/hari/Visiting_Research/Data/database/CPG_Island.csv', low_memory = False)
mafDataframe = pd.read_csv('/home/hari/Visiting_Research/Data/database/MAF_CSV.csv', low_memory = False)

#This list stores the methylation interval for all the entires in the Methylation database
methylationIntervalList = []

#traverse through the methylation database to find the methylation interval
for rowNumber in range(len(cpgDataframe)-1):
    #find the methylation sites
    methylationSites = (str)(cpgDataframe.loc[rowNumber, 'CpG_sites_hg19'])
    methylationSitesList = methylationSites.split(';')
    #consider the first and last site for interval calculation
    firstSite = methylationSitesList[0].split(":")
    lastSite = methylationSitesList[-1].split(":")
    chromosomeNumber = firstSite[0]
    firstSite = (int)(firstSite[1])
    lastSite = (int)(lastSite[1])
    #calculate the interval
    intervalStart = math.floor(firstSite - ((lastSite - firstSite) * 0.5))
    intervalEnd = math.ceil(lastSite + ((lastSite - firstSite) * 0.5))
    if intervalStart == intervalEnd:
        intervalStart-= 1
        intervalEnd+= 1
    #Append the methylation interval for the current sample in the database
    methylationIntervalList.append((str)((chromosomeNumber)+"_"+(str)(intervalStart)+"_"+(str)(intervalEnd)))

#Delete the last sample since it has only NA values
cpgDataframe.drop(index = 81037, inplace = True)

#Create a new column called Methylation_interval and assign the list as values
cpgDataframe['Methylation_interval'] = methylationIntervalList

#Key-value pairs where key is the chromosome number and the value is a list of indices in the Methylation database containing the samples of that chromosome.
#This is done to reduce the processing when we search for the methylation interval for a sample in the MAF database since we now only have to look at the samples of 
#that particular chromosomes in the methylation databse instead of the entire database until we find one. 
chromosomeMap = {(str)(chromosome) : [] for chromosome in range(1,23)}
chromosomeMap["X"] = []
chromosomeMap["Y"] = []
for rowNumber in range(len(cpgDataframe)):
    methylationInterval = cpgDataframe.loc[rowNumber, 'Methylation_interval']
    methylationIntervalList = methylationInterval.split('_')
    #chromosome is the key
    chromosome = methylationIntervalList[0]
    #Append the index to the list
    chromosomeMap[chromosome].append(rowNumber)

#This list will contain the list of all mutations from the MAF database lying in that corresponding methylation interval  
mutationsList = [ [] for i in range(len(cpgDataframe))]

#Traverse through each sample in the MAF database
for rowNumber in range(500):
    #Find the chromosome
    chromosome = (str)(mafDataframe.loc[rowNumber, 'Chromosome'])
    #Find the cell line
    cellLine = mafDataframe.loc[rowNumber, 'Tumor_Sample_Barcode']
    #if there are no samples in with that chromosome number or the cell line in the Methylation database, append an empty value and go to the next sample
    if ((cellLine not in cpgDataframe.columns.tolist()) or chromosome not in chromosomeMap.keys()):
        continue
    #Find the start and the end position of the sequence in the MAF database
    mafStart = mafDataframe.loc[rowNumber, 'Start_position']
    mafEnd = mafDataframe.loc[rowNumber, 'End_position']
    mafStart = (int)(mafStart)
    mafEnd = (int)(mafEnd)
    #Using the table Key-Value pairs we created earlier, we now only the samples with that chromosome number in the Methylation database 
    for index in chromosomeMap.get(chromosome):
        methylationInteval = cpgDataframe.loc[index, 'Methylation_interval']
        methylationIntervalList = methylationInteval.split('_')
        methylationStart = (int)(methylationIntervalList[1]) 
        methylationEnd = (int)(methylationIntervalList[2])
        #If either the start or the end is in the Methylation interval, then the MAF region overlaps with the methylation interval
        if((mafStart in range(methylationStart, methylationEnd+1)) or (mafEnd in range(methylationStart, methylationEnd+1))):
            #Find the methylation value corresponding to the cell line and append that to the list
            mutationsList[index].append(rowNumber);


#Write the output to a file. 
fileHandle = open('testRunResult.txt', 'a')
for i in range(len(cpgDataframe)):
    indexList = mutationsList[i]
    indexList = [str(indexList[j]) for j in range(len(indexList))]
    fileHandle.write(str(i)+"   ")
    fileHandle.write(','.join(indexList))
    fileHandle.write('\n')

fileHandle.close()
