// RLaB diary file: convert_laser_ADCcount_to_distance.diary.r. Opened Fri Nov 30 18:49:43 2018

HOME=811;
STEP1=2373;
STEP2=657;
STEP3=2372;
STEP4=657;
# not averaged, single readings
diff([HOME, STEP1])
        1562  
diff([STEP2, STEP1])
        1716  
diff([STEP2, STEP3])
        1715  
diff([STEP4, STEP3])
        1715  
9.378/1715
0.00546822157  
# 9.378 mm / ADC count diff =
diary()
