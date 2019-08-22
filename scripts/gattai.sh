FILE1=$1
FILE2=$2
NEW_FILE=$3

for i in $(seq 6 6 768); do
    head -n $i $FILE1 | tail -n 6 >> $NEW_FILE
    head -n $i $FILE2 |	tail -n	6 >> $NEW_FILE
done
