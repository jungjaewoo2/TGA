(awk '{print $1*4e-5, $2/4096*3.3 }' < new_amp_log_sweep_10Hz-3kHz_300mV_3s.dat   |tail  -n +2 |head -n -2 ;echo;\
 awk '{print $1*4e-5, $5/4096*3.3 }' < new_amp_log_sweep_10Hz-3kHz_300mV_3s.dat   |tail  -n +2 |head -n -2 ;echo;\
 awk '{print $1*4e-5, ($2 - $4)/4096*3.3 }' < new_amp_log_sweep_10Hz-3kHz_300mV_3s.dat   |tail  -n +2 |head -n -2 ;echo;\
 awk '{print $1*4e-5, $4/4096*3.3 }' < new_amp_log_sweep_10Hz-3kHz_300mV_3s.dat   |tail  -n +2 |head -n -2) |ap
