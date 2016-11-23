INPUT="$1"

while read line; do

   awk -v l="$line" '{if(substr(l,1,1) == ">"){if(substr(l,2,15) == $2){print ">"$3","$1","substr(l,2); exit}}else{print l; exit}}' "$INPUT"

done
