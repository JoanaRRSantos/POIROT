import os
import re
import pandas as pd
from Bio import Entrez
Entrez.email = "fetchmysequencefromgenbank@gmail.com"

class ProcessData:
    
    """
    This class will get and process the data given by the user in order to later be used by
    two other scripts depending on the need of the user. 
    The input file must be an excel with the name of the species that the user wants to search (for now)
    
    It is perfered if the excel has a column name saying Specie or Species.
    If NOT, the list of species MUST be in the first column of the excel.
    
    """

    def __init__(self, directory="", folderDirectoryToSaveFiles="", fileName="", dataset = "", headers ="", index =0, geneNamesToSearch=[], firstList=[], listOfUniqueSpecies=[], save=False, numberOfSeqs=20, IDs={}, min_length=100, max_length=2000):
        
        self.directory = directory #dictory of the folder where the initial input excel should be saved
        self.folderDirectoryToSaveFiles = folderDirectoryToSaveFiles #directory of the folder that files should be saved
        self.fileName = fileName #name of the excel with extension 
        self.dicExcel = {"Species":0, "Specie":1, "species":2, "specie":3, "declared species":4, "declared specie":5,"Declared species":6, "Declared specie":7, "Declared Species":8, "Declared Specie":9, "declared Species":10, "declared Specie":11, "Especies":12, "Espécies":13, "Especie":14, "Espécie":15, "especies":16, "espécies":17, "especie":18, "espécie":19, "especies declaradas":20, "especie declarada":21, "espécies declaradas":22, "espécie declarada":23, "Especies declaradas":24, "Especie declarada":25, "Espécies declaradas":26, "Espécie declarada":27, "Especies Declaradas":28, "Especie Declarada":29, "Espécies Declaradas":30, "Espécie Declarada":31, "especies Declaradas":32, "especie Declarada":33, "espécies Declaradas":33, "espécie Declarada":34}
        self.dataset = dataset #where the dataset is going to be stored after processing
        self.index = index #integer to know the column where the names of species is located
        self.headers = headers #if the excel given has a header, it will be stored here to add to the output excel
        self.geneNamesToSearch = geneNamesToSearch #list of genes the person wants to test in the search against the species given
        self.firstList = firstList #list of species with repited names 
        self.listOfUniqueSpecies = listOfUniqueSpecies
        self.save = save 
        self.yesDict = {"Yes":0, "yes":1, "sim":2, "Sim":3}
        self.numberOfSeqs = numberOfSeqs
        self.IDs= IDs
        self.min_length = min_length
        self.max_length = max_length

    def saveFilesFromUser(self):  
        """
        Fuction to know if the user wants to save the sequences found in his computer or not.
        If the user wants to, they should respond with a "yes" or "sim".
        """
        save_response = input("Do you want to save the sequences files in your PC?  (when done click ENTER)")
        if save_response in self.yesDict.keys():
            self.save=True

    def getDirectoryFromUser(self):
        """
        Function gets directory where the user has the input file and where they want to save the output files.
        The user must copy the directory of the folder(s) and paste it and hit ENTER.
        """
        if self.directory == "":
            folderDirectory = input("Tell me which directory you have the input excel: (when done click ENTER)")
            if os.path.exists(folderDirectory):
                self.directory = folderDirectory
            else:
                self.directory = os.getcwd()
            os.chdir(self.directory)
            if self.save:
                folderDirectoryToSaveFiles = input("Tell me which directory you want to save your files: (when done click ENTER)")
                if os.path.exists(folderDirectoryToSaveFiles):
                    self.folderDirectoryToSaveFiles = folderDirectoryToSaveFiles
                else:
                    self.folderDirectoryToSaveFiles = self.directory                                                                                                                                                    

    def getFileNameFromUser(self):
        """
        Function gets the name of an excel where the user saved all the species names that they want to search.
        The user must write the name of the excel without the termination. 
        """
        fileName = input("Tell me the name of the excel which has the species names:  (when done click ENTER)")
        fileName += ".xlsx"
        try:
            if os.path.exists(fileName):
                self.fileName = fileName
        except NameError:
            print("The input file needs to be saved in the same folder as the directory given where the excel is saved")

    def getGenesFromUser(self):
        """
        Function gets the names of the genes that the user wants to search about 
        The user must write the name and click ENTER
        When there is no more genes, the user must simply click ENTER for the program to start
        """
        while True:
            line = input("What are the names of the genes you want to search: (when done click ENTER)")
            if line:
                self.geneNamesToSearch.append(line)
            else:
                break   
            
    def getNumberOfSequencesFromUser(self):
        """
        Function to know how many sequences the user wants to save from the records
        The user should respond here with any positive number, if no number is given, 20 sequences is set by default
        """
        number_response = input("How many sequences do you want to save from the records?")
        try:
            self.numberOfSeqs = int(number_response) if number_response.strip() else 20
        except ValueError:
            self.numberOfSeqs = 20

    def getLengthOfSeqsFromUser(self):
        """
        Function to know the length of the sequences retrived from the record
        The user should respond here with any positive number, if no number is given, the minimum is set to 100bp and the maximum is set to 2000bp
        """
        min_length = input("Enter the minimum length of the sequence (e.g., 100): ")
        max_length = input("Enter the maximum length of the sequence (e.g., 2000): ")
        try:
            self.min_length = int(min_length) if min_length.strip() else 100
        except ValueError:
            self.min_length = 100
        try:
            self.max_length = int(max_length) if max_length.strip() else 2000
        except ValueError:
            self.max_length = 2000


    def processExcel(self):
        """
        Function gets directory and excel name given by the user and processes to get a dataset with the organization needed for the program to work.
        If the excel given has an header, the column with the species can be anywhere.
        If the excel does not have an header, the column with the species MUST be in the first column.
        """
        file = os.path.join(self.directory, self.fileName)
        df = pd.read_excel(file) #da maneira que esta deduz que vai ter header  # can also index sheet by name or fetch all sheets
        df = df.T.reset_index().T.reset_index(drop=True)
        self.dataset = df
        self.headers = df.iloc[0].tolist()
        for i in range(len(self.headers)):
            if self.headers[i] in self.dicExcel.keys():
                self.index = i
                self.dataset = df.iloc[1:]    
    def GetIDsFromGenBank(self):
        """
        Function that searches for the genes for every species given and will save in an dictionary of dicgtonaries the information of {species : {gene: [list of ids found]}}
        This information will later be used by either fetchMySequences.py or ____.py 
        """
        self.firstList = self.dataset[self.index].values.tolist()
        self.firstList = [item.rstrip() if item.endswith(" ") else item for item in self.firstList]
        [self.listOfUniqueSpecies.append(species) for species in self.firstList if species not in self.listOfUniqueSpecies]
        for species in self.listOfUniqueSpecies:
            self.IDs[species] = {}
            for gene in self.geneNamesToSearch:
                tempSearch = species + "[Orgn] AND " + gene
                handle = Entrez.esearch(db="nucleotide", term= tempSearch, idtype="acc", retmax=self.numberOfSeqs)
                record = Entrez.read(handle)
                self.IDs[species][gene] = record["IdList"]

        

if __name__ == "__main__":
    teste = ProcessData()
    teste.saveFilesFromUser()
    teste.getDirectoryFromUser()
    teste.getFileNameFromUser()
    teste.getGenesFromUser()
    teste.getNumberOfSequencesFromUser()
    teste.processExcel()
    teste.GetIDsFromGenBank()
