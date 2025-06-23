# set -x
# Manual calibration of coarse motion tester

FILE=/home/pi/rlab/lib.r/specs.r
LCORR=$(grep LCORR $FILE | tail -1 | cut -d ';' -f 1 | cut -d = -f 2 | tr -d ' ')

# MIN=-15;
# MIN=`echo $MIN|cut -d '.' -f 1`
CENTER=0
RANGE=15
MIN=$(($CENTER - $RANGE))
MAX=$(($CENTER + $RANGE))
CORR=$(
	/usr/bin/zenity \
		--width=400 \
		--title "Course motion tester calibration" \
		--scale \
		--value=$CENTER \
		--min-value=$MIN \
		--max-value=$MAX \
		--text "Move the slider to change calibration" \
		2>/dev/null
)

if [ "$?" = "0" ]; then
	OP=$(cat /var/www/html/Operator.txt)
	DATE=$(date)
	LCORR=$(($LCORR + $CORR))

	echo "LCORR = $LCORR;# modified by $OP on $DATE" >>$FILE
fi
