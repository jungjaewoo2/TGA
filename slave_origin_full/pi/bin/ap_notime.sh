(awk '{print NR*4e-5, $1/4096*3.3 - 1.65 }' < $1  |tail  -n +2 |head -n -2 ;echo;\
 awk '{print NR*4e-5, $2/4096*3.3 - 1.65 }' < $1  |tail  -n +2 |head -n -2 ;echo;\
 awk '{print NR*4e-5, $3/4096*3.3 - 1.65 }' < $1  |tail  -n +2 |head -n -2 ;echo;\
 awk '{print NR*4e-5, $4/4096*3.3 - 1.65 }' < $1  |tail  -n +2 |head -n -2 ;echo;\
 awk '{print NR*4e-5, $5/4096*3.3 }' < $1  |tail  -n +2 |head -n -2) |ap -w $2
