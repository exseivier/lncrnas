#!/usr/bin/env bash

BATCH_SIZE=$1
FILE=$2

#	DISCLAIMERS
[ -f $FILE ] || echo "[ERROR!] - File $FILE was not found"
[ -f $FILE ] || exit 1
echo "$FILE was found"
FILE_SIZE=$(cat $FILE | wc -l)
FILE_HEADERS_NUMBER=$(grep ">" $FILE | wc -l)
[ $FILE_HEADERS_NUMBER -eq $(echo "$FILE_SIZE / 2" | bc) ] || echo "[ERROR!] - Lines in file $FILE do not match with header numbers"
[ $FILE_HEADERS_NUMBER -eq $(echo "$FILE_SIZE / 2" | bc) ] || exit 1
echo "Testing sequences file size and header numbers... [OK]"
[ $(echo "$BATCH_SIZE % 2" | bc) -eq 0 ] || echo "[ERROR!] - Batch size is not even"
[ $(echo "$BATCH_SIZE % 2" | bc) -eq 0 ] || exit 1
echo "Testing if batch size is even... [OK]"
echo ""

PIVOT=$FILE_SIZE
PID_NUM=$RANDOM
> idonotknowwhy.txt
> stats.dgs.${PID_NUM}.txt
for ((i=1; i <= $FILE_SIZE; i = i + $BATCH_SIZE));
do
	echo "tail -n $PIVOT $FILE | head -n $BATCH_SIZE > tmp.${PID_NUM}"
	if [ $(echo "$PIVOT - $BATCH_SIZE" | bc -l) -le 0 ];
	then
		BATCH_SIZE=$PIVOT
	else
		PIVOT=$(echo "$PIVOT - $BATCH_SIZE" | bc -l)
	fi
	echo "tail -n $PIVOT $FILE | head -n $BATCH_SIZE > tmp.${PID_NUM}"
	echo "cat tmp.${PID_NUM}"
	
	echo "RNAfold -i tmp.${PID_NUM} -o idonotknowwhy.txt --noPS"
	echo "cat *.fold >> stats.dgs.${PID_NUM}.txt"
	echo "rm *.fold"
done
rm idonotknowwhy.txt
