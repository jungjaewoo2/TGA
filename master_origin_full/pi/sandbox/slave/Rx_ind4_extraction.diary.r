// RLaB diary file: Rx_ind4_extraction.diary.r. Opened Sun Apr 28 17:26:11 2019

k1=Rdc_extract_start[1]
       23105  
k4=Rdc_extract_start[2]
       39432  
kstep=int((k4-k1)/3)
        5442  
ind4=k1+(0:4)*kstep
       23105         28547         33989         39431         44873  
plot([[ind4[1]:ind4[5]]',d[ind4[1]:ind4[5];]])
