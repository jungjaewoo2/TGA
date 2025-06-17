// RLaB diary file: ADC_coarse_motion_laser_stepper_conersion.diary.r. Opened Fri Nov 30 17:59:55 2018

# graphically read from ap plot
# laser RS232 results 8.579 mm, 9.378 mm, 9.378 mm, 9.378 mm
# for STEP1 STEP2 STEP3 STEP4 travel
#
# actual direct laser (RS232) readings
# -006.733 +001.786
# -007.574
# +001.786
# -007.573
# the above are laser readings
# home - step1
-006.733 -001.786
      -8.519  
# home - step1 = -8.519 
diff([-006.733, +001.786])
       8.519  
# step1 - step2
diff([-007.574, +001.786])
        9.36  
# step3 - step2
diff([-007.574, +001.786])
        9.36  
# step3 - step4
diff([-007.573, +001.786])
       9.359  
# CORRECTION: direct laser readings (not from plot) are:
# 8.519, 9.36 , 9.36  9.359 

# stepper step size:
9.36/740*2*1e3
  25.2972973  
9.36/740*2*1e3
  25.2972973  
9.359/740*2*1e3
  25.2945946  
8.519/672.5*2*1e3
   25.335316  
# versus IBM's 25.4166 
diary();
