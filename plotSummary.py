import pandas as pd
import matplotlib.pyplot as plt

mafDataFrame = pd.read_csv('/home/hari/Visiting_Research/Data/database/MAF_CSV.csv', low_memory=False);
cpgDataFrame = pd.read_csv('/home/hari/Visiting_Research/Data/database/CPG_Island.csv', low_memory=False);

cancerCellLineMap = {}
seenCellLines = []

for rowNumber in range(len(mafDataFrame)):
    cellLine = (str)(mafDataFrame.loc[rowNumber, 'Tumor_Sample_Barcode'])
    if cellLine not in (cpgDataFrame.columns.tolist()):
        continue
    if cellLine in seenCellLines:
        continue
    seenCellLines.append(cellLine)
    cellLineList = cellLine.split("_")
    cancerType = "_".join(cellLineList[1:])
    if cancerType in cancerCellLineMap:
        cancerCellLineMap[cancerType]+= 1
    else:
        cancerCellLineMap[cancerType] = 1

count = 0

labels = []
values = []

for key,val in cancerCellLineMap.items():
    labels.append(key)
    values.append(val)

print(labels)
print(values)
figureObject, axesObject = plt.subplots()
figureObject.set_facecolor('black')

axesObject.pie(values, labels = labels, autopct='%1.2f', startangle=90, radius=20)

axesObject.axis('equal')

plt.show()