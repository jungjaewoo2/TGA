// RLaB diary file: debugging_inverted_curves.diary.r. Opened Sun Mar 10 09:24:53 2019

"stepper_twang_rev3.for_debugging.r"
stepper_twang_rev3.for_debugging.r"using prod_Mar9_inverted_5.dat "
using prod_Mar9_inverted_5.dat diff(ind)
 matrix columns 1 thru 9
          28          2693            28          2694            28          2693            28          2693            28  

 matrix columns 10 thru 18
        2693            27          2444            28          2693            28          2693            28          2693  

 matrix columns 19 thru 27
          28          1362          5440             3          5439             3          5439             3          5439  

# but BEMF_extract_start = ind[find((diff(ind) >2500) && (diff(ind) <3000 ))];
# so 2444 is ignored
#
# 2. BEMF_extract_start count both positive and negative pulses
#    It is better to select only positive pulses like:
# define new ind3
ind3=find(ddiff2);
find((diff(ind3) >24) && (diff(ind3) <32 )) 
           1             3             5             7             9  
ind3[max(find((diff(ind3) >24) && (diff(ind3) <32 )))]
       29328  

# verify time on plot for "prod_Mar9_inverted_5.dat"
ind3[max(find((diff(ind3) >24) && (diff(ind3) <32 )))]*4e-5
     1.17312  
# or previous one
ind3[max(find((diff(ind3) >24) && (diff(ind3) <32 ))) -1 ]*4e-5
     0.95656  
(BEMF_extract_start)*4e-5
     0.20464       0.31348       0.42236        0.5312       0.64004       0.84772       0.95656        1.0654  
# the same time stamps so skip value shouls stay the same
diary();
