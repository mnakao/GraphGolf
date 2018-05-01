for t in 1.865207
	 do
for i in $(seq 1 20)
do
      ./GraphGolf -g 2 -f ../data/n16d5.random.edges -n 1000000 -s $i -w $t -c $t  > log.$t.$i.txt &
done
wait
done
