# ./lcdi2c -b 1 -l -x 0 -y 0  $(date)
# ./lcdi2c -i -b 1 -l -x 6 -y 0  $(date)

 ./lcdi2c -i 
while true; do 
	./lcdi2c -b 1 -l -x 0 -y 0  $(date)

	sleep 1; 
done
