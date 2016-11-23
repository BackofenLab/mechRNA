#Adds last column from a file to the end of another file
#Usage: add_log_fold_change.sh [SOURCE_FILE] [FILE] [INT_INDEX1] [INT_INDEX2]
source="$1"
file="$2"
index1="$3"
index2="$4"

while read g; do
  key=$(echo "$g" | awk -v i="$index2" '{print $i}')

  awk -F'[\t]' -v i="$index1" -v key="$key" -v line="$g" '{if($i == key){print line"\t"$NF}}' "$source"

done < "$file"
