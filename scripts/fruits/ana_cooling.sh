#!/bin/bash

TEXT=result.txt

for C in 1 10 100 1000 10000; do
  for i in $(seq 1 10); do
    grep "ASPL gap" log.C$C.$i.txt | awk '{print $16}'
done | sort > $C
done

paste 1 10 100 1000 10000 > $TEXT
echo "" >> $TEXT

for C in 1 10 100 1000 10000; do
  for i in $(seq 1 10); do
    grep Steps log.C$C.$i.txt | awk '{print $2}'
  done | sort -n > $C
done

paste 1 10 100 1000 10000 >> $TEXT

rm -f 1 10 100 1000 10000


