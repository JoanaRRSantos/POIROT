import os
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "fetchmysequencefromgenbank@gmail.com"

tempSearch=["Rhodiola internal transcribed spacer"]
for term in tempSearch:
    print(term)
    handle = Entrez.esearch(db="nucleotide", term=term, idtype="acc", retmax=300)
    record = Entrez.read(handle)
    print("im getting IDs")
    with open("outFile", 'w') as w:
        for ID in record["IdList"]:
            print(ID)
            results = Entrez.efetch(db="nucleotide", id=ID, rettype="fasta")
            fetch_record = SeqIO.read(results, "fasta")
            results.close()
            SeqIO.write(fetch_record, "current_seq.fasta", "fasta")
            for line in open('current_seq.fasta'):
                w.write(line)
            os.remove("current_seq.fasta")