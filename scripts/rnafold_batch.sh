#!/usr/bin/env bash

BATCH_SIZE=$1
FILE=$2

#	DISCLAIMERS
[ -f $FILE ] || echo "[ERROR!] - File $FILE was not found"
[ -f $FILE ] || exit 1
echo ""
echo "$FILE was found... [OK]"
FILE_SIZE=$(cat $FILE | wc -l)
FILE_HEADERS_NUMBER=$(grep ">" $FILE | wc -l)
[ $FILE_HEADERS_NUMBER -eq $(echo "$FILE_SIZE / 2" | bc) ] || echo "[ERROR!] - Lines in file $FILE do not match with header numbers"
[ $FILE_HEADERS_NUMBER -eq $(echo "$FILE_SIZE / 2" | bc) ] || exit 1
echo "Testing sequences file size and header numbers... [OK]"
[ $(echo "$BATCH_SIZE % 2" | bc) -eq 0 ] || echo "[ERROR!] - Batch size is not even"
[ $(echo "$BATCH_SIZE % 2" | bc) -eq 0 ] || exit 1
echo "Testing if batch size is even... [OK]"
echo "Processing RNA fold predictions and calculating stability energy..."

PIVOT=$FILE_SIZE
PID_NUM=$RANDOM
mkdir fold.${PID_NUM}
#cd fold.${PID_NUM}
> fold.${PID_NUM}/stats.dgs.${PID_NUM}.txt
for ((i=1; i <= $FILE_SIZE; i = i + $BATCH_SIZE));
do
	if [ $(echo "$PIVOT - $BATCH_SIZE" | bc -l) -le 0 ];
	then
		BATCH_SIZE=$PIVOT
	else
		PIVOT=$(echo "$PIVOT - $BATCH_SIZE" | bc -l)
	fi
	tail -n $PIVOT $FILE | head -n $BATCH_SIZE > fold.${PID_NUM}/tmp.${PID_NUM}
	#cat tmp.${PID_NUM}
	RNAfold -i fold.${PID_NUM}/tmp.${PID_NUM} --noPS >> fold.${PID_NUM}/stats.dgs.${PID_NUM}.txt
	#cat fold.${PID_NUM}/*.fold >> fold.${PID_NUM}/stats.dgs.${PID_NUM}.txt
	#rm fold.${PID_NUM}/*.fold
done
