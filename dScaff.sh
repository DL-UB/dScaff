#!/bin/bash

# Input files: 
# -a (assembly) = draft assembly in FASTA format
# -q (query) = gene.fna file with gene sequences from reference genome or ranked_queries.fasta output of ranked queries SubSequencesExtractor.sh script
# -d (dataset) = ncbi_dataset.tsv containing all genes in reference genome or coordinates_dataset.csv output of ranked queries SubSequencesExtractor.sh script
# -gq for gene queries strategy
# -rq for ranked queries strategy



# Define the help method
get_help() {
  echo " "
  echo "Usage: $(basename $0) [-h] [-a FILE] [-q FILE] [-d FILE] [OPTION]"
  echo " "
  echo "Options:"
  echo "  -h, --help			Display this help message"
  echo "  -a, --assembly		Draft assembly in FASTA format"
  echo "  -q, --query			gene.fna file with gene sequences from reference genome or ranked_queries.fasta output of ranked queries SubSequencesExtractor.sh script"
  echo "  -d, --dataset			ncbi_dataset.tsv containing all genes in reference genome or coordinates_dataset.csv output of ranked queries SubSequencesExtractor.sh script"
  echo "  -gq, --gene_queries		Perform gene queries strategy"
  echo "  -rq, --ranked_queries		Perform ranked queries strategy"
  echo " "
  echo "Use files from working directory or provide absolute path to input files."
  echo " "
  exit 0
}

# Check if no arguments provided
if [ $# -eq 0 ]; then
  echo " "
  echo "Error: No arguments provided."
  echo " "
  get_help
fi

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  get_help
fi

# Initialize variables
assembly=""
query=""
dataset=""
strategy=""

while [ ! -z "$1" ]; do
  case "$1" in
     --assembly|-a)
         shift
         echo " "
         echo "Input draft assembly is: $1"
         assembly=$1
         if [ -z "$1" ]; then
            echo " "
            echo "Error: Missing argument for draft assembly file."
            echo " "
            exit 1
         fi         
         ;;
     --query|-q)
         shift
         echo " "
         echo "Input reference query sequences is: $1"
         query=$1
         if [ -z "$1" ]; then
            echo " "
            echo "Error: Missing argument for reference query sequences."
            echo " "
            exit 1
         fi  
         ;;
     --dataset|-d)
        shift
        echo " "
        echo "Input dataset is: $1"
        dataset=$1 
        if [ -z "$1" ]; then
            echo " "
            echo "Error: Missing argument for dataset."
            echo " "
            exit 1
         fi         
        ;;
     --gene_queries|-gq)
         echo " "
         echo "Gene queries strategy selected."
         strategy="gene_queries"
         ;;
     --ranked_queries|-rq)
         echo " "
         echo "Ranked queries strategy selected."
         strategy="ranked_queries"
         ;;
     *)
        get_help
        ;;
  esac
  shift
done

# Validate required inputs
if [ -z "$assembly" ]; then
    echo " "
    echo "Error: Missing input draft assembly file"
    echo " "
    exit 1
fi
if [ -z "$query" ]; then
    echo " "
    echo "Error: Missing input reference query file"
    echo " "
    exit 1
fi
if [ -z "$dataset" ]; then
    echo " "
    echo "Error: Missing input dataset file"
    echo " "
    exit 1
fi

# Check if a strategy was selected
if [ -z "$strategy" ]; then
    echo " "
    echo "Error: No strategy selected. Please use -gq for gene queries or -rq for ranked queries."
    echo " "
    exit 1
fi



START=$(date +%s)


### spinner function
spinner() {
    local pid=$1
    local delay=0.1
    local spinstr='|/-\'
    while [ "$(ps a | awk '{print $1}' | grep $pid)" ]; do
        local temp=${spinstr#?}
        printf "  %c  " "$spinstr"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b\b"
    done
    printf "    \b\b\b\b"
}








###
# First part --> filter genes by chromosomes
###
CWD="$(pwd)"

echo " "
echo "Identifing chromosomes ..."
echo " "
awk -F "\t" 'NR>1{print $4}' $dataset | sort | uniq > chromosomes.txt

sed -i '/^$/d' chromosomes.txt

subdir="$(basename -- ${assembly%%.*})"
mkdir $subdir

echo "Selecting queries for each chromosome ..."
echo " "
cat chromosomes.txt | while read line 
do 
awk -v line=$line -F "\t" '{if ($4 == line) print $0}' $dataset > ./$subdir/$line.tsv
done

mv chromosomes.txt $subdir


echo "Filtering datasets ..."
echo " "
cp query_filtering.R $subdir
cp parameters_dScaff.txt $subdir
cd $subdir
Rscript query_filtering.R
rm parameters_dScaff.txt


###
# Second part --> blastn between each query of interest and draft assembly
###

cd $CWD
#lungimi contiguri
awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($1,2,100) ","; } $0 !~ ">" {c+=length($0);} END { print c; }' $assembly > lungimi_contiguri.csv

echo "Preparing BLAST database ..."
mkdir assembly_database
cp $assembly assembly_database
cd assembly_database
makeblastdb -in *.fasta -dbtype nucl
cd $CWD/$subdir/

echo "Performing BLAST for selected queries ..."
#echo " "
ls -d ./*/ | while read line
do 

#[[ ! -d "line" ]] && continue
 
cd $line
echo " "
echo " "
echo "Entered " $line


awk -F "," 'NR > 1 {print $2":"$3"-"$4}' *_filtered.csv > query_headers.txt
sed -i 's/"//g' query_headers.txt

cd $CWD/$subdir
####################################################################################################################################
grep ">" $CWD/$query | sed 's/>//g' > headers_list_db.txt
####################################################################################################################################
cd $line
cat query_headers.txt | while read subquery; do grep $subquery ../headers_list_db.txt >> query_headers_db.txt; done

####################################################################################################################################
seqtk subseq $CWD/$query query_headers_db.txt > queries_of_interest.fasta
####################################################################################################################################
mkdir queries

awk -F ">| " '/^>/ {s=$2".query.fasta"}; {print > s}' queries_of_interest.fasta

mv *.query.fasta ./queries

cd ./queries
file_count=$(ls *.fasta | wc -l)
current_file=0
cd ..

threads=$(nproc --all)

for i in queries/*.query.fasta
do

((current_file++))
percent=$(( 100 * current_file / file_count ))

blastn -db $CWD/assembly_database/*.fasta -query $i -out ${i%.query.fasta}".lucru.csv" -outfmt 6 -num_threads $threads

# build progress bar
    progress_bar=""
    for ((j=0; j<percent; j+=2)); do
        progress_bar="${progress_bar}#"
    done
    for ((j=percent; j<100; j+=2)); do
        progress_bar="${progress_bar}-"
    done
# print progress bar
echo -ne "[${progress_bar}] ${percent} %\r"

done

for i in queries/*.query.fasta; do lung=$( echo $i | sed -e '1d' $i | wc -c ); length=$( expr $lung - 1 ); awk -v Length=$length -F "\t" '{ FS = OFS = "\t" } {print $0,Length}' ${i%.query.fasta}".lucru.csv" > ${i%.query.fasta}".csv"; done


rm -r queries/*.query.fasta queries/*.lucru.csv
#rm -r queries

cd $CWD/$subdir
rm headers_list_db.txt

done



###
# Third part --> map draft assembly contigs using the genes of interest
###
echo " "
echo " "
echo "Indexing and mapping contigs ..."

ls -d ./*/ | while read line 
do 

#[[ ! -d "line" ]] && continue

cd $CWD/$subdir/$line

#echo "Entered " $line

mv *_distances_filtered.csv query_filtered.csv

cd queries
find . -name '*' -size 0 -print0 | xargs -0 rm 2> /dev/null


done

cd $CWD
cp contigs_mapping.R $subdir
cp parameters_dScaff.txt $subdir
mv lungimi_contiguri.csv $subdir 
cd $subdir

Rscript contigs_mapping.R 2> /dev/null &
spinner $!
#echo "Done!"
rm parameters_dScaff.txt
rm lungimi_contiguri.csv

############################
cd $CWD
grep ">" $assembly > headers_assembly.txt
############################
cd $subdir

ls -d ./*/ | while read line 
do
cd $line

mkdir tmp
mv *.csv tmp
mv *.txt tmp
mv *.fasta tmp

#######################################
	if [ -d "chromosome" ]; then
	cd chromosome/
	if [ $strategy == "gene_queries" ]; then
	rm *orientation.csv
	fi
##################################################################################################################
	sed -i 's/"//g' *selected_contigs.csv 
	sed -i 's/,/ /g' *selected_contigs.csv
#Pentru fiecare fisier care se termina cu selected_contigs.csv, il parcurg cu while, apoi extrag headerele de interes in baza fisierului cu all headers.
	for i in *selected_contigs.csv; do cat $i | while read contig index; do grep $contig $CWD/headers_assembly.txt >> ${i%_selected_contigs.csv}".selected.headers.txt"; done; done

#Adaug o coloana care contine cromozomul si scaffoldul.
	for i in *selected_contigs.csv; do awk -v loc="${i%_selected_contigs.csv}" '{print $0, loc}' $i >> contigs_chr_loc.txt; done

#Elimin ">".
	sed -i 's/>//g' *selected.headers.txt

#Extrag contigurile, dar nu sunt in ordine.
	for i in *.selected.headers.txt; do seqtk subseq $CWD/$assembly $i > ${i%.selected.headers.txt}".new1.assembly.fasta"; done

#Sparg fisierul fasta ce contine contiguri.
	awk -F ">| " '/^>/ {s=$2".fna"}; {print > s}' *.new1.assembly.fasta

#Am ordonat contigurile
	cat contigs_chr_loc.txt | while read contig loc chr; do for i in *.fna; do if [[ ${i%.fna} == $contig ]]; then cat $i >> $chr".dScaff.assembly1.fasta"; fi; done; done
	rm -r *.fna
	rm -r *.new1.assembly.fasta

	for chr in *.dScaff.assembly1.fasta
	do
	awk '/^>/{gsub(/^>/,">"i++" ");}1' i=1 $chr > ${chr%.dScaff.assembly1.fasta}"_dScaff_assembly.fasta"
	sed -i "s/>/>${chr%.dScaff.assembly1.fasta}"_"/g" ${chr%.dScaff.assembly1.fasta}"_dScaff_assembly.fasta"
	done
	rm -r contigs_chr_loc.txt
	rm -r *.dScaff.assembly1.fasta
	rm -r *.selected.headers.txt
	######################deni 31 oct
	rm -r *_selected_contigs.csv
	######################deni 31 oct
	sed -i 's/_chromosome//g' *_dScaff_assembly.fasta
##################################################################################################################

	fi
#######################################
#######################################
		if [ -d "scaffolds" ]; then
		cd scaffolds 2> /dev/null
		ls -d ./*/ | while read scaff 1> /dev/null
		do
		cd $scaff
		if [ $strategy == "gene_queries" ]; then
		rm *orientation.csv
		fi 

##################################################################################################################
		if test -f *selected_contigs.csv; then

		sed -i 's/"//g' *selected_contigs.csv 
		sed -i 's/,/ /g' *selected_contigs.csv
#Pentru fiecare fisier care se termina cu selected_contigs.csv, il parcurg cu while, apoi extrag headerele de interes in baza fisierului cu all headers.
		for i in *selected_contigs.csv; do cat $i | while read contig index; do grep $contig $CWD/headers_assembly.txt >> ${i%_selected_contigs.csv}".selected.headers.txt"; done; done

#Adaug o coloana care contine cromozomul si scaffoldul.
		for i in *selected_contigs.csv; do awk -v loc="${i%_selected_contigs.csv}" '{print $0, loc}' $i >> contigs_chr_loc.txt; done

#Elimin ">".
		sed -i 's/>//g' *selected.headers.txt

#Extrag contigurile, dar nu sunt in ordine.
		for i in *.selected.headers.txt; do seqtk subseq $CWD/$assembly $i > ${i%.selected.headers.txt}".new1.assembly.fasta"; done

#Sparg fisierul fasta ce contine contiguri.
		awk -F ">| " '/^>/ {s=$2".fna"}; {print > s}' *.new1.assembly.fasta

#Am ordonat contigurile
		cat contigs_chr_loc.txt | while read contig loc chr; do for i in *.fna; do if [[ ${i%.fna} == $contig ]]; then cat $i >> $chr".dScaff.assembly1.fasta"; fi; done; done
		rm -r *.fna
		rm -r *.new1.assembly.fasta

		for chr in *.dScaff.assembly1.fasta
		do
		awk '/^>/{gsub(/^>/,">"i++" ");}1' i=1 $chr > ${chr%.dScaff.assembly1.fasta}"_dScaff_assembly.fasta"
		sed -i "s/>/>${chr%.dScaff.assembly1.fasta}"_"/g" ${chr%.dScaff.assembly1.fasta}"_dScaff_assembly.fasta"
		done
		rm -r contigs_chr_loc.txt
		rm -r *.dScaff.assembly1.fasta
		rm -r *.selected.headers.txt
		######################deni 31 oct
		rm -r *_selected_contigs.csv
		######################deni 31 oct
		sed -i 's/_scaffold//g' *_dScaff_assembly.fasta

##################################################################################################################
		cp *.fasta ../
		fi
		cd - 1> /dev/null
		done
		
		line_name=$( echo $line | cut -c3- )
		cat *.fasta > ${line_name%/}"_dScaff_assembly.fa"
		rm *.fasta
		cat *.fa > ${line_name%/}"_dScaff_assembly.fasta"
		rm *.fa
		fi
#######################################
cd $CWD/$subdir

find . -type f -name *_dScaff_assembly.fasta -exec cat > dScaff.fasta {} +

done

rm query_filtering.R
rm contigs_mapping.R
rm chromosomes.txt
rm *.tsv

cd $CWD
#mv assembly_database $subdir
rm -r assembly_database
rm headers_assembly.txt

find . -name queries -type d -exec rm -rf {} +


echo " "
echo "Finished !"
### End
echo " "
END=$(date +%s)
seconds=$(( $END - $START ))
minutes=$(echo "scale=1; $seconds / 60" | bc)
hours=$(echo "scale=1; $seconds / 3600" | bc)

if [[ $seconds -ge 3600 ]] 
then
echo "dScaff ran for $hours hours (that's $minutes minutes or $seconds seconds)."
elif [[ $seconds -ge 60 ]]
then
echo "dScaff ran for $minutes minutes (that's $seconds seconds)."
else
echo "dScaff ran for $seconds seconds."
fi

echo " "

# CWD este dScaff folder
# subdir este folder cu numele asamblarii
# line este folder cu numele chromozomului

