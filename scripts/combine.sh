if [ $# -ne 3 ]; then
    echo $0 "FILE1 FILE2 (nodes*based_degree/2/g)"
    exit 0
fi

FILE1=$1
FILE2=$2
START=$3  # (nodes*based_degree/2/g)
STEP=$START
END=$(wc -l $FILE1 | awk '{print $1}')
rm -f $NEW_FILE
for i in $(seq $START $STEP $END); do
    head -n $i $FILE1 | tail -n $STEP
    head -n $i $FILE2 | tail -n $STEP
done
