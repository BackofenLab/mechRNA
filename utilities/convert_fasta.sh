INPUT="$1"
OUTPUT="$2"
LEN="$3"

awk -v len="$LEN" 'BEGIN{header = ""; seq = ""}{if($0 ~ ">"){if(header != "" && length(seq) >= len){print header; print seq}; header = $0; seq = ""}else{seq = seq$0}}END{print header; print seq}' "$INPUT"  | awk -F'[>: ]' '{if($0 ~ ">" && (length($6) == 1 || length($6) == 2)){print ">"substr($2,1,15)";"substr($11,1,15)";"$17; getline; print $0}}' |paste -d '$' - - |sort -t";" -k 2|tr '$' '\n' > "$OUTPUT.$LEN"
