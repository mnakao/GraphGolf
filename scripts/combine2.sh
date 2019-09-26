if [ $# -ne 3 ]; then
    echo $0 "FILE1 FILE2 (based_nodes*based_degree/2)"
    exit 0
fi

FILE1=$1
FILE2=$2
START=$3
STEP=$START
END=$(wc -l $FILE1 | awk '{print $1}')
for i in $(seq $START $STEP $END); do
    head -n $i $FILE1 | tail -n $STEP
    head -n $i $FILE2 | tail -n $STEP
done
