#!/bin/zsh

########################################################################
###### Data processing pipeline for ITS2 Metabarcoding
###### by Alexander Keller (University of Wuerzburg)
###### for INSIGNIA, www.insignia-bee.eu (contact: Alice Pinto)
########################################################################

# define variables
seqfilter=$2
ps=$3
vsearch=$4
threads=$5

cd $1
classificationOnly=($(grep "classificationOnly" config.txt | cut -f2 -d"="))

  #extracting files
  find . -name '*.gz' -print0 | xargs -0 -I {} -P $threads gunzip {}

  mkdir -p logs

  for f in *_R1_*.fastq; do

    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
    s=$(cut -d_ -f1 <<< "$f")
    p=$(cut -d_ -f2 <<< "$f")
  	total=$(grep "^@M02" $f | wc -l)


    echo " "
    echo "===================================="
    echo "Processing sample $s"
    echo "===================================="
    $vsearch --fastq_mergepairs  $f --reverse $r --fastq_minovlen 20 --fastq_maxdiffs 10 --fastqout $s.merge.fq --fastq_eeout --relabel R1+2-$s- --threads $threads  2> logs/vsearch.m.$s.log

    var="$(grep "Merged" logs/vsearch.m.$s.log)"
    echo "$s : $var" | tee -a logs/_merging.log

    $vsearch --fastq_filter $s.merge.fq \
          --fastq_maxee 1 \
          --fastq_minlen 200 \
          --fastq_maxlen 550 \
          --fastq_maxns 0 \
          --fastaout $s.mergefiltered.fa \
          --fasta_width 0 --threads $threads 2> logs/vsearch.mf.$s.log

    var="$(grep "sequences kept" logs/vsearch.mf.$s.log)"
    echo "$s : $var" | tee -a logs/_filter.log

    $vsearch --fastq_truncee 1.5 --fastq_filter $f --fastq_minlen 200 --fastaout $s.trunc.fa --relabel R1-$s- --threads $threads 2> logs/vsearch.tf.$s.log

    var="$(grep "sequences kept" logs/vsearch.tf.$s.log)"
    echo "$s : $var" | tee -a logs/_truncfilter.log

  done

  echo " "
  echo "===================================="
  echo "ASV generation and mapping"
  echo "===================================="


  cat *mergefiltered.fa > all.merge.fasta
  cat *trunc.fa > all.trunc.fasta

  # cleanup
  mkdir -p raw
  mkdir -p tmp

  mv *.fastq raw/
  mv *.fq tmp/
  mv *.fa tmp/

  #tar

  $vsearch --derep_fulllength all.merge.fasta \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.merge.derep.uc \
    --output all.merge.derep.fa --threads $threads 2> logs/_derep.log

  $vsearch --cluster_unoise  all.merge.derep.fa --sizein --sizeout --centroids zotus.merge_chim.fa --threads $threads 2> logs/_unoise.log

  grep "Clusters" logs/_unoise.log
  grep "Singletons" logs/_unoise.log

  $vsearch --sortbysize zotus.merge_chim.fa --output zotus.merge_sorted.fa --threads $threads 2>  logs/_sort.log

  $vsearch --uchime3_denovo zotus.merge_sorted.fa --abskew 16 --nonchimeras zotus.merge.fa --threads $threads 2>  logs/_uchime.log
  cp zotus.merge.fa  zotus.merge.fa.bak
  $seqfilter zotus.merge.fa.bak  --ids-rename='s/^.*$/'$1'-$COUNT/' --out zotus.merge.fa
  grep "Found" logs/_uchime.log

  # [TODO] make decision whether to take merged or R1 here at some point

  ### create community table
  cat all.merge.fasta |  sed "s/^>R1+2-\([a-zA-Z0-9-]*\)\-\([0-9]*\)/>R1+2-\1_\2;barcodelabel=\1;/g" > all.merge.bc.fasta
  cat all.trunc.fasta |  sed "s/^>R1-\([a-zA-Z0-9-]*\)\-\([0-9]*\)/>R1-\1_\2;barcodelabel=\1;/g" > all.trunc.bc.fasta

  $vsearch --usearch_global all.merge.bc.fasta --db zotus.merge.fa --strand plus --id 0.97 --uc map.merge.uc --threads $threads 2> logs/_mapping.log

  grep "Matching" logs/_mapping.log

  #python2.7 $p/uc2otutab.py map.trunc.uc > zotu_table.trunc.txt
  python2.7 $ps/uc2otutab.py map.merge.uc > asv_table.merge.txt

##### create taxonomy

  echo " "
  echo "===================================="
  echo "Taxonomic classification"
  echo "===================================="

# $vsearch
refDBs=($(grep "refdb" config.txt | cut -f2 -d"=" | sed 's/\"//g'))
hieDBs=($(grep "hiedb" config.txt | cut -f2 -d"=" | sed 's/\"//g'))

threshold=97

echo ",kingdom,phylum,order,family,genus,species" > taxonomy.vsearch
#echo ",kingdom,phylum,order,family,genus,species" > taxonomy.blast

countdb=0
cp  zotus.merge.fa zotus.direct.$countdb.uc.nohit.fasta
prevDB=$countdb

for db in "${refDBs[@]}"
  do :
    countdb=$((countdb+1))
    echo "\n\n#### Direct VSEARCH Classification level: $countdb";
    $vsearch --usearch_global zotus.direct.$prevDB.uc.nohit.fasta --db $db --id 0.$threshold --uc zotus.direct.$countdb.uc --fastapairs zotus.direct.$countdb.fasta --strand both --threads $threads 2>  logs/_direct.$countdb.log

    grep "^N[[:space:]]" zotus.direct.$countdb.uc | cut -f 9 > zotus.direct.$countdb.uc.nohit
    $seqfilter zotus.merge.fa --ids zotus.direct.$countdb.uc.nohit --out zotus.direct.$countdb.uc.nohit.fasta
    cut -f 9,10 zotus.direct.$countdb.uc  | grep -v "*" | sed "s/[A-Za-z0-9]*;tax=//" >> taxonomy.vsearch
    prevDB=$countdb
  done

echo "\n\n#### Hierarchical VSEARCH classification";

$vsearch --sintax zotus.direct.$countdb.uc.nohit.fasta -db $hieDBs -tabbedout zotus.uc.merge.nohit.sintax -strand plus -sintax_cutoff 0.9 -threads $threads 2>  logs/_sintax.log

cut -f1,4 zotus.uc.merge.nohit.sintax | sed -E -e "s/\_[0-9]+//g" -e "s/,s:.*$//"  >> taxonomy.vsearch

#BLAST local DBs
#countdb=0
#cp  zotus.merge.fa zotus.blast.$countdb.uc.nohit.fasta
#prevDB=$countdb
#touch zotus.blast.hits
#touch taxonomy.blast

#for db in "${refDBs[@]}"
 # do :
  #  countdb=$((countdb+1))
   # echo "\n\n#### Direct BLAST Classification level: $countdb";
    #makeblastdb -in $db -parse_seqids -blastdb_version 5 -dbtype nucl
  #  blastn  -outfmt '6 qseqid sseqid length pident qcovs' -max_target_seqs 1  -query  zotus.blast.$prevDB.uc.nohit.fasta -subject $db -perc_identity $threshold -qcov_hsp_perc 90 -num_threads $threads > zotus.blast.$countdb.out
  #  cut -f1 -d"	" zotus.blast.$countdb.out | cut -f1 >> zotus.blast.hits
  #  $seqfilter zotus.merge.fa --ids-exclude --ids  zotus.blast.hits --out zotus.blast.$countdb.uc.nohit.fasta
  #  prevDB=$countdb
  #  cut -f 1,2 zotus.blast.$countdb.out >> taxonomy.blast

  #done

#BLAST web nt instead of sintax here?
#Alex
#sed -i .bak -e "s/c:.*,o:/o:/g" -e "s/[A-Za-z0-9]*;tax=//" -e "s/	/,/" taxonomy.vsearch
#Rufino
sed -i.bak -e "s/c:.*,o:/o:/g" -e "s/[A-Za-z0-9]*;tax=//" -e "s/	/,/" taxonomy.vsearch
#Alex
#sed -i .bak -e "s/c:.*,o:/o:/g" -e "s/[A-Za-z0-9]*;tax=//" -e "s/	/,/" taxonomy.blast
#Rufino
#sed -i.bak -e "s/c:.*,o:/o:/g" -e "s/[A-Za-z0-9]*;tax=//" -e "s/	/,/" taxonomy.blast
