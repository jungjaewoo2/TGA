// UPDATED on Tue Jul 16 11:41:05 KST 2019
// DEBUG2=1; // Testing at home without lasers 
DEBUG2=0; 
UNDER = 20;
OVER = -5;
#RCAL=0.63;
#RCAL=0.29;
RCAL=0.0;
##################################
# @27.25 degC
RCAL1=0.21;
RCAL2=0.0;
RCAL1_20C=0.17;
RCAL2_20C=0.0;
##################################
RDCMAX=15.9;
RDCMIN=13.9;
TAUmin=2.3;
TAUmax=5.9;
FMAX=100;
FMIN=70;
RESRISEmin = 12.0;
RESRISEmax = 19.3;
# TADDR1="28-02131ace25aa"; 
# TADDR1="28-02131acde2aa";
# TADDR1="28-011927954fa2"; # work station No.2
###########    TADDR1="28-02131a8169aa"; #  work station No.3
# TADDR1="28-01142fc8aefe"; # Home
pBEMF =[ 0.0113279432,    63.4832683,    109.919378,   -2.36119721,  0.00120285811 ]';
pLASER=[0.0103179959,      64.68188,    56.4212137,    -1.4772254,  -0.00291975192]';
pLASER=[0.0109,      64.68188,    65,    -1.7,  -0.0029]';
#pLASER=[ 0.0290637238,    82.3755559,    32.3201767,   -1.71585538,  0.00336292393 ]';
# temporary for inverted magnets
# pLASER=[0.0103179959,      64.68188,    56.4212137,    -1.4772254,  -0.00291975192]';
# pLASER=[ 0.0101240054,    89.3421555,    33.8407838,  -0.587874212,  0.0024213385 ]';
#pLASER=[ 0.010452566,    80.2664382,    38.7644328,   -0.80269388,  0.00124523479 ]';
# for TGA rev3 [azg] Tue Jul 16 11:40:36 KST 2019
pLASER=[0.00930804502,     68.647364,    88.3570708,   -1.18213193,   -3.05968487 ]';
pLASER=[0.00930804502,     68.647364,    137,   -1.5,   0.0 ]';
pLASER= [0.00940908385,     72.949328,    83.2823856,   0.473120464,  -0.000443014547 ]';
TRAVELmin=1200;
GAINmin=2.45;
SKEWGAINmin=0.92;
SKEW_TRAVELmin=193; # NOT CLEAR IF REQ'D
SKEW_ANGLEmin=1.1 ;
# LCORR = -4.43180556;
# LCORR = -19.49; # tester #2
LCORR = 60.11; # tester #3
LCORR = 48; # tester #3
LCORR = 49;# modified by Andrew on Tue Aug 13 04:29:49 KST 2024
LCORR = 48;# modified by Andrew on Tue Aug 13 04:31:40 KST 2024
LCORR = 47;# modified by Andrew on Tue Aug 13 07:53:20 KST 2024
LCORR = 49;# modified by Andrew on Tue Aug 13 07:55:29 KST 2024
LCORR = 44;# modified by SBK_3 on Sat Dec 24 15:00:48 KST 2022
LCORR = 41;# modified by SBK_3 on Sat Dec 24 15:09:10 KST 2022
LCORR = 47;# modified by SBK_3 on Sat Dec 24 15:16:00 KST 2022
LCORR = 48;# modified by SBK_3 on Sat Dec 24 15:18:39 KST 2022
LCORR = 46;# modified by  on Mon Sep  9 10:15:31 KST 2024
LCORR = 48;# modified by SBK_3 on Mon Sep  9 10:30:49 KST 2024
