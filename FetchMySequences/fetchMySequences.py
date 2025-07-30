
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 14:02:28 2023

@author: Joana Santos
"""
from Main import ProcessData
import os
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "fetchmysequencefromgenbank@gmail.com"

class fetchMySequences(ProcessData):
    
    """
    This class will fetch from the processed data the IDs obtain and process them as the user wants.
    It can download to diferent folders the gene sequences for each species and analise the results (long search) or
    It can determine if there is sequences of a certain gene for a certain species and the number of sequences (quick search).
    
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def SearchSequencesInGenBank(self):
        """
        Function that searches for the genes for every species given and will give as an output a new excel with the results added as well as save the sequences found if the user wants it to
        """
        if self.save:
            main_output_filename="Results_"+self.fileName+".xlsx"
            for gene in self.geneNamesToSearch:
                tempDict={}    
                for species in self.listOfUniqueSpecies:
                    tempDirectory = os.path.join(self.folderDirectoryToSaveFiles, species)
                    if not os.path.exists(tempDirectory):
                        os.mkdir(tempDirectory)
                    os.chdir(tempDirectory)
                    output_filename = species + "_" + gene + ".fasta"
                    with open(output_filename, 'a') as w:
                        if self.IDs[species][gene] != []:
                            tempDict[species]="yes"
                        for ID in self.IDs[species][gene]:
                            results = Entrez.efetch(db="nucleotide", id=ID, rettype="fasta", retmode="text")
                            try:
                                fetch_record = SeqIO.read(results, "fasta")
                            except ValueError as e:
                                print("Error reading record:", e)
                                continue  
                            finally:
                                    results.close()
                            SeqIO.write(fetch_record, "current_seq.fasta", "fasta")
                            for line in open('current_seq.fasta'):
                                w.write(line)
                            os.remove("current_seq.fasta")
                        ColumnToAdd = ["no"] * len(self.firstList)
                        ColumnToAdd2 = [0] * len(self.firstList)
                        for species in range(len(self.firstList)):
                            for key in tempDict.keys():
                                if self.firstList[species] == key:
                                    ColumnToAdd[species] = "yes"
                                    ColumnToAdd2[species] = len(self.IDs[self.firstList[species]][gene])
                self.dataset[gene] = ColumnToAdd
                self.dataset["number of sequences of " + gene] = ColumnToAdd2
            self.dataset.columns.values[0:len(self.headers)] = self.headers 
            os.chdir(self.folderDirectoryToSaveFiles)
            self.dataset.to_excel(main_output_filename)
            print("Your results where saved succesfully, please check your folder for the new excel")
        else:
            output_filename2="Results_"+self.fileName+".xlsx"
            for gene in self.geneNamesToSearch:
                tempDict = {}
                for species in self.listOfUniqueSpecies:
                    if self.IDs[species][gene] != []:
                        tempDict[species]="yes"
                ColumnToAdd = ["no"] * len(self.firstList)
                ColumnToAdd2 = [0] * len(self.firstList)
                for species in range(len(self.firstList)):
                    for key in tempDict.keys():
                        if self.firstList[species] == key:
                            ColumnToAdd[species] = "yes"
                            ColumnToAdd2[species] = len(self.IDs[self.firstList[species]][gene])
                self.dataset[gene] = ColumnToAdd
                self.dataset["number of sequences of " + gene] = ColumnToAdd2
            self.dataset.columns.values[0:len(self.headers)] = self.headers 
            self.dataset.to_excel(output_filename2)
            print("Your results where saved succesfully, please check your folder for the new excel")

if __name__ == "__main__":
    teste = fetchMySequences()
    teste.saveFilesFromUser()
    teste.getDirectoryFromUser()
    teste.getFileNameFromUser()
    teste.getGenesFromUser()
    teste.getNumberOfSequencesFromUser()
    teste.processExcel()
    teste.GetIDsFromGenBank()
    teste.SearchSequencesInGenBank()
