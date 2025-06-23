Rdc_20degC = function ( Rdc, temp){
// function returns resistance value at 20 deg C
//
// Rdc -> uncorrected resistance at temp
// RCAL -> wire resistance at 20 deg C
// temp -> temperature 
// Mon Jul 15 18:27:45 PDT 2019
//////////////////////
alpha20=40e-4;

Rdc_20deg = Rdc /(1+ alpha20*(temp-20)); 


	return Rdc_20deg ;
};

