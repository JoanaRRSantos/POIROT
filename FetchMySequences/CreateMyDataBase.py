from Main import ProcessData
import os
from Bio import Entrez, SeqIO
import time
import xml.etree.ElementTree as ET
from bitstring import BitArray  # Make sure you're using this library for binary conversion


class createDataBase(ProcessData):
    
    """
    This class will fetch from the processed data the IDs obtain and process them as the user wants.
    It can download to diferent folders the gene sequences for each species and analise the results (long search) or
    It can determine if there is sequences of a certain gene for a certain species and the number of sequences (quick search).
    
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.desired_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        self.rank_abbr = ['k', 'p', 'c', 'o', 'f', 'g', 's']
        
    def SearchSequencesAndTaxonomyInGenBank(self):
        max_retries = 2   
        for gene in self.geneNamesToSearch:
            filename = "Database_" + gene + ".txt"
            sequences_dict = {}
            for species in self.listOfUniqueSpecies:
                print("Searching species: ", species)
                if self.IDs[species][gene] != []:
                    tempSearch = f"{species}[Organism]"
                    for attempt in range(max_retries):
                        try:
                            handle = Entrez.esearch(db="taxonomy", term=tempSearch, retmax=1)
                            record = Entrez.read(handle)
                            handle.close()
                            taxonomy_handle = Entrez.efetch(db="taxonomy", id=record["IdList"][0], retmode="xml")
                            taxonomy_record = Entrez.read(taxonomy_handle)
                            taxonomy_handle.close()
                            taxonomy = [item["ScientificName"] for item in taxonomy_record[0]["LineageEx"] if item["Rank"] in self.desired_ranks]
                            #if taxonomy does not have species in ranks
                            if len(taxonomy) == 6:
                                tempSpeciesWanted = str(species).replace(" ", "_")
                                taxonomy.append(tempSpeciesWanted)
                            #if taxonomy has species in ranks
                            elif len(taxonomy) == 7:
                                tempSpeciesWanted = str(taxonomy[-1]).replace(" ", "_")
                                taxonomy[-1]=tempSpeciesWanted
                            tempTaxWanted = ','.join(f"{abbr}:{name}" for abbr, name in zip(self.rank_abbr, taxonomy))
                            break
                        except Exception as e:
                            print(f"Error fetching taxonomy for {species} (attempt {attempt + 1}): {e}")
                            time.sleep(2)
                    else:
                        print(f"Failed to retrieve taxonomy for {species} after {max_retries} attempts")
                        continue
                        
                    failed_ids = set(self.IDs[species][gene])
                    for attempt in range(max_retries):
                        to_retry = failed_ids.copy()
                        failed_ids.clear()
                        for ID in to_retry:
                            print(ID, species)
                            try:
                                results = Entrez.efetch(db="nucleotide", id=ID, rettype="gb", retmode="text")
                                gb_data = SeqIO.read(results, "genbank")
                                organism = gb_data.annotations['organism']
                                temp_species = ' '.join(organism.split()[:2]).replace(" ","_")
                                record_taxonomy = gb_data.annotations['taxonomy']
                                record_taxonomy.append(temp_species)
                                results.close()                        
                                temp_taxonomy_elements = [item.split(":")[1] for item in tempTaxWanted.split(",")]
                                if all(taxon in record_taxonomy for taxon in temp_taxonomy_elements):
                                    if ID not in sequences_dict:
                                        sequences_dict[ID] = {"taxonomy": tempTaxWanted, "sequence": str(gb_data.seq)}  
                            except Exception as e:
                                failed_ids.add(ID)
                                print(f"Error fetching sequence for ID {ID} (attempt {attempt + 1}): {e}")
                                time.sleep(2)
                        if not failed_ids:
                            break
                    if failed_ids:
                        print(f"Failed to retrieve sequences for these IDs after {max_retries} attempts: {failed_ids}")
        
            with open(filename, "a") as f:
                for ID, data in sequences_dict.items():
                    if self.min_length <= len(data["sequence"]) <= self.max_length:
                        tempID=str(ID).split('.')[0]
                        taxonomy = data["taxonomy"]
                        sequence = data["sequence"]
                        f.write(f">{tempID};tax={taxonomy}\n{sequence}\n")
                    else:
                        print(f"Sequence for ID {ID} removed because it is not in the range of {self.min_length} to {self.max_length} base pairs.")
            print("The database of " + gene + " has been created")
        print("All databases are available in your input file folder")    
                        
                    
if __name__ == "__main__":
    teste = createDataBase()
    teste.getDirectoryFromUser()
    teste.getFileNameFromUser()
    teste.getGenesFromUser()
    teste.getNumberOfSequencesFromUser()
    teste.getLengthOfSeqsFromUser()
    teste.processExcel()
    teste.GetIDsFromGenBank()
    teste.SearchSequencesAndTaxonomyInGenBank()
                        
