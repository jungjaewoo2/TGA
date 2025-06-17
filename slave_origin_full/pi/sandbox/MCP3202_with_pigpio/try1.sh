(awk '{print $1*4e-5, $2/4096*3.3 }' < fromRPi2_1st_try.dat |tail  -n +2 |head -n -2 ;echo;\
 awk '{print $1*4e-5, $3/4096*3.3 }' < fromRPi2_1st_try.dat |tail  -n +2 |head -n -2 ;echo;\
 awk '{print $1*4e-5, $5/4096*3.3 }' < fromRPi2_1st_try.dat |tail  -n +2 |head -n -2 ;echo;\
 awk '{print $1*4e-5, $4/4096*3.3 }' < fromRPi2_1st_try.dat |tail  -n +2 |head -n -2) |ap
