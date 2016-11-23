# Fixes result file in cases when 4 sub-optimals do not exist
# Usage: fill_missing.sh [FILE]
result="$1"
good=0
refseq2=""

read line
while [ $good ]
do
   refseq=$(echo "$line" | awk -F '[\t]' '{print $1}')

   if [ "$refseq" == "" ]; then
           good=1
           break
   fi

   echo "$line"

   for i in {1..4}
   do
	read line
	refseq2=$(echo "$line" | awk -F '[\t]' '{print $1}')

	if [ "$refseq2" == "" ]; then
	   good=1
           break
        fi

	if [ "$refseq" != "$refseq2" ]; then
	   for (( j=$i; j<=4; j++ ))
	   do
   	      echo "$refseq"$'\t'NA$'\t'NA$'\t'NA$'\t'NA$'\t'NA
	   done
	   break
        else
	   echo "$line"
	fi
   done
   if [ "$refseq" == "$refseq2" ]; then
	read line
   fi
done < "$result"
