(awk '{print $1*4e-5, $2/4096*3.3 }' < new_amp_10Hz_square_vol_0.25.dat  |tail  -n +2 |head -n -2 ;echo;\
 awk '{print $1*4e-5, $5/4096*3.3 }' < new_amp_10Hz_square_vol_0.25.dat  |tail  -n +2 |head -n -2 ;echo;\
 awk '{print $1*4e-5, ($2 - $4)/4096*3.3 }' < new_amp_10Hz_square_vol_0.25.dat  |tail  -n +2 |head -n -2 ;echo;\
 awk '{print $1*4e-5, $4/4096*3.3 }' < new_amp_10Hz_square_vol_0.25.dat  |tail  -n +2 |head -n -2) |ap
