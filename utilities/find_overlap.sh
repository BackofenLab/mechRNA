#Checks if an IntaRNA hit overlaps the GraphProt Peak
#Usage: find_overlap.sh [INT_DIST] [FILE]

dist=$1
file=$2

awk -F'\t' -v d=$dist '{if($2 < (151+d) && $3 > (151-d)){print $0}}' "$file"
