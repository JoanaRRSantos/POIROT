# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 15:38:10 2023

@author: Joana Santos
"""

file1 = open('speciesPDFeEXCELcomDUPLICADOScomAMOSTRAS.txt', 'r')
Lines = file1.readlines()
SpeciesDict = {}
for line in Lines:
    line = line.strip()
    if line not in SpeciesDict.keys():
        SpeciesDict[line] = 1
    else: 
        SpeciesDict[line] += 1

#print(SpeciesDict)

count=0
for values in SpeciesDict.values():
   count += values
print("Number of species incialy in file: " ,count)

count=0
with open ("ListOfSpeciesForDatabase.txt", "w", encoding="utf-8") as f:
    for species in SpeciesDict.keys():
        count+=1
        f.write(species)
        f.write('\n')
f.close()
print("Number of species without repeats in file: " ,count)


