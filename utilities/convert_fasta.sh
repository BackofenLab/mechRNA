INPUT="$1"
OUTPUT="$2"
LEN="$3"

awk -v len="$LEN" 'BEGIN{header = ""; seq = ""}{if($0 ~ ">"){if(header != "" && length(seq) >= len){print header; print seq}; header = $0; seq = ""}else{seq = seq$0}}END{print header; print seq}' "$INPUT"  | awk -F'[>: ]' '{if($0 ~ ">" && (length($7) == 1 || length($7) == 2)){print $0; getline; print $0}}' > "$OUTPUT.$LEN"
