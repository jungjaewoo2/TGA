// RLaB diary file: skew_polarity_and_sequence_test.diary.r. Opened Thu Jun  6 18:07:28 2019

skew_pol=sign(diff(Ldist)[3]/diff(Current)[3])
           1  
# skew must be much(?) bigger in the last two intervals
abs(diff(Ldist)[3]) > abs(diff(Ldist)[1])
           1  
if( abs(diff(Ldist)[3]) > abs(diff(Ldist)[1]) ){ SKEW_SEQ="Normal";} else SKEW_SEQ="Reversed"
if( abs(diff(Ldist)[3]) > abs(diff(Ldist)[1]) ){ SKEW_SEQ="Normal";} 
if( abs(diff(Ldist)[3]) < abs(diff(Ldist)[1]) ){ SKEW_SEQ="Reversed";} 
if( skew_pol >0 ) { SKEW_POL="OK";}
if( skew_pol <0 ) { SKEW_POL="Reversed";}
SKEW_SEQ
NormalSKEW_SEQ
NormalSKEW_POL
OKdiary();
