# use $1 as an identifier
# STATUS BROKEN: scp fails (empty file)
# USAGE: $0 7
#
sudo ~/bin/prodMCP3202 > /tmp/prod"$1"_slave.dat &
ssh tga.local  'sudo /home/pi/bin/prodMCP3202 > /tmp/prod'$1'_tga.dat &' &
ssh gen.local  '/usr/bin/play /home/pi/bin/5Hz_1ms_TGA.flac vol 0.3' 
sleep 7
# scp tga.local:/tmp/prod${1}* /tmp
# scp slave.local:/tmp/prod${1}* /tmp
scp /tmp/prod"$1"_slave.dat  tga.local:/tmp/prod_slave_remote.dat
