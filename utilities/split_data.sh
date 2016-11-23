#Splits input data into separate files
#Usage: split_data.sh [INT_NUM_FILES] [FILE]

num_files=$1
file=$2

total=$(wc -l $file | awk '{print $1}')

num_lines=$(($total/$num_files))
num_lines=$(($num_lines+1))

split -d -l $num_lines $file input
mv input00 input0
mv input01 input1
mv input02 input2
mv input03 input3
mv input04 input4
mv input05 input5
mv input06 input6
mv input07 input7
mv input08 input8
mv input09 input9

