import pandas as pd 
import numpy as np 
import math
from intervalTree import Interval
from intervalTree import IntervalTree
import sys


cpgDataframe = pd.read_csv('/home/hari/Documents/Visiting_Research/Data/database/CPG_Island.csv', low_memory = False)
mafDataframe = pd.read_csv('/home/hari/Documents/Visiting_Research/Data/database/MAF_CSV.csv', low_memory = False)


methylationIntervalList = []

chromosomeMap = {(str)(chromosome) : None for chromosome in range(1,23)}
chromosomeMap["X"] = None
chromosomeMap["Y"] = None

testNode = IntervalTree()
print(sys.getrecursionlimit())

for rowNumber in range(len(cpgDataframe)-1):
    methylationSites = (str)(cpgDataframe.loc[rowNumber, 'CpG_sites_hg19'])
    methylationSitesList = methylationSites.split(';')
    firstSite = methylationSitesList[0].split(":")
    lastSite = methylationSitesList[-1].split(":")
    chromosomeNumber = firstSite[0]
    firstSite = (int)(firstSite[1])
    lastSite = (int)(lastSite[1])
    intervalStart = math.floor(firstSite - ((lastSite - firstSite) * 0.5))
    intervalEnd = math.ceil(lastSite + ((lastSite - firstSite) * 0.5))
    curInteval = Interval(intervalStart, intervalEnd, rowNumber)
    chromosomeMap[chromosomeNumber] = testNode.insertNode(chromosomeMap[chromosomeNumber], curInteval)
    


cpgDataframe.drop(index = 81037, inplace = True)

methylationValuesList = []


for rowNumber in range(len(mafDataframe)):
    print(rowNumber)
    chromosome = (str)(mafDataframe.loc[rowNumber, 'Chromosome'])
    cellLine = mafDataframe.loc[rowNumber, 'Tumor_Sample_Barcode']
    methylationValuesIndices = []
    if ((cellLine not in cpgDataframe.columns.tolist()) or chromosome not in chromosomeMap.keys()):
         methylationValuesList.append(methylationValuesIndices)
         continue
    mafStart = mafDataframe.loc[rowNumber, 'Start_position']
    mafEnd = mafDataframe.loc[rowNumber, 'End_position']
    mafStart = (int)(mafStart)
    mafEnd = (int)(mafEnd) 
    testNode.intervalSearch(chromosomeMap[chromosome], Interval(mafStart, mafEnd, -1), methylationValuesIndices)
    methylationValuesList.append(methylationValuesIndices)

