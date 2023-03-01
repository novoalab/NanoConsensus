################################################
### Script to extract full-length reads' IDs ###
################################################

#Read inputs from the command line: 
bed=$1
chr=$2
start=$3
end=$4
output_name=$5

#Parse bed file with awk and output results:
echo $chr
awk -v chr="$chr" -v start="$start" -v end="$end" '{if($1==chr && $2<=start && $3>=end) {print $4}}' $bed > $output_name.$chr.read_IDs.txt 
