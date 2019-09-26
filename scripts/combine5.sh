if [ $# -ne 6 ]; then
    echo $0 "FILE1 FILE2 FILE3 FILE4 FILE5 (based_nodes*based_degree/2)"
    exit 0
fi

FILE1=$1
FILE2=$2
FILE3=$3
FILE4=$4
FILE5=$5
START=$6
STEP=$START
END=$(wc -l $FILE1 | awk '{print $1}')
for i in $(seq $START $STEP $END); do
    head -n $i $FILE1 | tail -n $STEP
    head -n $i $FILE2 | tail -n $STEP
    head -n $i $FILE3 | tail -n $STEP
    head -n $i $FILE4 | tail -n $STEP
    head -n $i $FILE5 | tail -n $STEP
done
