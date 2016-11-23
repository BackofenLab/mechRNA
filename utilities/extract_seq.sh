#Extract sequences from fasta file given a list of IDs (gene [tab] id;id;id...)
source=$1

while read id
do
awk -v F="$id" '{if($1 ~ F){printf "%s", RS $0};  RS=">"; FS="\n"}' "$source"
done
