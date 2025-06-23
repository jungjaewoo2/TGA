#!/bin/bash
DEBUG=0
if [ $DEBUG = "1" ]; then
	echo "DEBUG MODE"
fi

DIALOG="`which zenity` --width 400"
TITLE="--title="
TEXT="--text="
ENTRY="--entry "
ENTRYTEXT="--entry-text "
# MENU="--list --print-column=1 --column=Pick --column=Info"
MENU="--list --hide-column=1 --print-column=1 --column=Pick --column=Info"
YESNO="--question "
MSGBOX="--info "
SCALE="--scale "
TITLETEXT="Date and Time Setting Tool"



while [ "$SETCHOICE" != "Exit" ]; do
DAY="`date +%d`"
MONTH="`date +%m`"
YEAR="`date +%Y`"
MINUTE="`date +%M`"
SETMINUTE=$MINUTE
HOUR="`date +%H`"
SETHOUR=$HOUR

SETCHOICE=`$DIALOG --height 300 $TITLE"$TITLETEXT" $MENU $TEXT"Set Year Date and Time\n\nTime=$HOUR:$MINUTE\nDate=$MONTH-$DAY-$YEAR\n\n" Exit "Quit" SETTIME "Set Current Time" SETDATE "Set Current Date" 2>/dev/null` 
SETCHOICE=`echo $SETCHOICE | cut -d "|" -f 1`

if [ "$SETCHOICE" = "SETTIME" ]; then

HOUR="`date +%H`"
HOUR=`echo $HOUR | sed -e 's/^0//g'`

SETHOUR=`$DIALOG $TITLE"$TITLETEXT" $SCALE --value=$HOUR --min-value=0 --max-value=23 $TEXT"Move the slider to the correct Hour" 2>/dev/null`
if [ "$?" = "0" ]; then

if [ "${#SETHOUR}" = "1" ]; then
SETHOUR="0$SETHOUR"
fi

MINUTE="`date +%M`"
MINUTE=`echo $MINUTE | sed -e 's/^0//g'`
fi

SETMINUTE=`$DIALOG $TITLE"$TITLETEXT" $SCALE --value=$MINUTE --min-value=0 --max-value=59 $TEXT"Move the slider to the correct Minute" 2>/dev/null`
if [ "$?" = "0" ]; then

if [ "${#SETMINUTE}" = "1" ]; then
SETMINUTE="0$SETMINUTE"
fi


if [ $DEBUG = "1" ]; then
echo  "date $MONTH$DAY$SETHOUR$SETMINUTE$YEAR"
echo "hwclock --systohc"
else
sudo date $MONTH$DAY$SETHOUR$SETMINUTE$YEAR
sudo hwclock --systohc
fi
fi
fi

if [ "$SETCHOICE" = "SETDATE" ]; then
DAY="`date +%d`"
DAY=`echo $DAY | sed -e 's/^0//g'`
MONTH="`date +%m`"
MONTH=`echo $MONTH | sed -e 's/^0//g'`
YEAR="`date +%Y`"

if CAL=`/usr/bin/zenity --calendar \
--title="Select a Date" \
--text="Click on a date to select that date." \
--day=${DAY} --month=${MONTH} --year=${YEAR} 2>/dev/null`
  then echo $CAL
MONTH=`echo ${CAL} | cut -d "/" -f 1`
if [ "${#MONTH}" = "1" ]; then
	MONTH="0$MONTH"
fi
DAY=`echo ${CAL} | cut -d "/" -f 2`
if [ "${#DAY}" = "1" ]; then
	DAY="0$DAY"
fi
YEAR=`echo ${CAL} | cut -d "/" -f 3`

  else echo "No date selected"
fi


MINUTE="`date +%M`"
HOUR="`date +%H`"

if [ $DEBUG == "1" ];then
echo "date $MONTH$DAY$HOUR$MINUTE$YEAR"
echo "hwclock --systohc"
else
sudo date $MONTH$DAY$HOUR$MINUTE$YEAR
sudo hwclock --systohc
fi
fi


done


exit 0
