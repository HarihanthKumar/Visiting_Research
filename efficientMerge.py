from os import sep
import pandas as pd 
import numpy as np 
import math
#The next two lines import Interval and Interval tree classes that will be used to store the methylation intervals in an efficient way 
from intervalTree import Interval
from intervalTree import IntervalTree
import sys

#load the databases as pandas dataframe
cpgDataframe = pd.read_csv('/home/hari/Visiting_Research/Data/database/CPG_Island.csv', low_memory = False)
mafDataframe = pd.read_csv('/home/hari/Visiting_Research/Data/database/MAF_CSV.csv', low_memory = False)

#create a dictionary that stores the top most root element for each chromosome and initialize each root value to null
chromosomeMap = {(str)(chromosome) : None for chromosome in range(1,23)}
chromosomeMap["X"] = None
chromosomeMap["Y"] = None

#create an instance of IntervalTree class to call the methods in the class
testNode = IntervalTree()

#print recursion stack size to check how deep can the recursion call be (python has an initial value of 1000)
print(sys.getrecursionlimit())

#Traverse each entry 
for rowNumber in range(len(cpgDataframe)-1):
    #find all the methylation sites and store them as a list.
    methylationSites = (str)(cpgDataframe.loc[rowNumber, 'CpG_sites_hg19'])
    methylationSitesList = methylationSites.split(';')
    #The methylation sites are already sorted in the ascending order of the sequence number and hence no sorting is needed.
    #Consider the first and the last site from the list.  
    firstSite = methylationSitesList[0].split(":")
    lastSite = methylationSitesList[-1].split(":")
    #Separate the chromosome number, first site and last site sequence numbers.  
    chromosomeNumber = firstSite[0]
    firstSite = (int)(firstSite[1])
    lastSite = (int)(lastSite[1])
    #Find the interval start and interval end
    intervalStart = math.floor(firstSite - ((lastSite - firstSite) * 0.5))
    intervalEnd = math.ceil(lastSite + ((lastSite - firstSite) * 0.5))
    if intervalStart == intervalEnd:
        intervalStart-= 1
        intervalEnd+= 1
    #create a new Interval object using the interval start, end and the index (position in the CpG database) values. 
    #This index value will be the one used later to access the methylation values
    curInteval = Interval(intervalStart, intervalEnd, rowNumber)
    #Add the interval to the self balancing interval tree and update the root interval. 
    chromosomeMap[chromosomeNumber] = testNode.insertNode(chromosomeMap[chromosomeNumber], curInteval)
    

#This is done to remove the final entry in the CpG database since it only has NA values.
cpgDataframe.drop(index = 81037, inplace = True)

#This list will contain the list of all mutations from the MAF database lying in that corresponding methylation interval  
mutationsList = [ [] for i in range(len(cpgDataframe))]

#Traverse each entry in the MAF database to find the list of all the overlapping indices from the CpG database
for rowNumber in range(len(mafDataframe)):
    #Find the chromosome number
    chromosome = (str)(mafDataframe.loc[rowNumber, 'Chromosome'])
    #Find the cell line of the entry
    cellLine = mafDataframe.loc[rowNumber, 'Tumor_Sample_Barcode']
    #If the cell line is not in the CpG database or if the chromosome is 'M' skip (Add an empty list).
    if ((cellLine not in cpgDataframe.columns.tolist()) or chromosome not in chromosomeMap.keys()):
        continue
    #find the starting and ending positions of the mutation 
    mafStart = mafDataframe.loc[rowNumber, 'Start_position']
    mafEnd = mafDataframe.loc[rowNumber, 'End_position']
    mafStart = (int)(mafStart)
    mafEnd = (int)(mafEnd) 
    #Search the interval tree and add the mutation entry to all the overalapping methylation intervals (The list
    # is passed to the intervalSearch function and the mutation entries are added there)
    testNode.intervalSearch(chromosomeMap[chromosome], Interval(mafStart, mafEnd, -1), mutationsList, rowNumber)


#Write the output to a file. 
fileHandle = open('testRunResult.txt', 'a')
for i in range(len(cpgDataframe)):
    indexList = mutationsList[i]
    indexList = [str(indexList[j]) for j in range(len(indexList))]
    fileHandle.write(str(i)+"   ")
    fileHandle.write(','.join(indexList))
    fileHandle.write('\n')

fileHandle.close()