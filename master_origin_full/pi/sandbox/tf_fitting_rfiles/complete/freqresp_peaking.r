// For a second order system, the damping ratio (zeta) determines the peaking in
//  the frequency domain (Peaking exists only for Zeta<0.707)
//
freqresp_peaking = function(zeta)
{
################################
# same formula as:
#  1/(2*zeta*sqrt(1 - zeta^2))
################################
	# allow calc negative "peaking"
	# if (zeta[1] > 0.707 ){ error ("Zeta<0.707 No peaking possible");}
return dB=-20*log10(2*zeta.*(1-zeta.^2).^0.5);
};
