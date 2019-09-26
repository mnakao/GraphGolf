if [ $# -ne 11 ]; then
    echo $0 "FILE[1-10] (based_nodes*based_degree/2)"
    exit 0
fi

FILE1=${1}
FILE2=${2}
FILE3=${3}
FILE4=${4}
FILE5=${5}
FILE6=${6}
FILE7=${7}
FILE8=${8}
FILE9=${9}
FILE10=${10}
START=${11}
STEP=$START
END=$(wc -l $FILE1 | awk '{print $1}')
for i in $(seq $START $STEP $END); do
    head -n $i $FILE1  | tail -n $STEP
    head -n $i $FILE2  | tail -n $STEP
    head -n $i $FILE3  | tail -n $STEP
    head -n $i $FILE4  | tail -n $STEP
    head -n $i $FILE5  | tail -n $STEP
    head -n $i $FILE6  | tail -n $STEP
    head -n $i $FILE7  | tail -n $STEP
    head -n $i $FILE8  | tail -n $STEP
    head -n $i $FILE9  | tail -n $STEP
    head -n $i $FILE10 | tail -n $STEP
done
