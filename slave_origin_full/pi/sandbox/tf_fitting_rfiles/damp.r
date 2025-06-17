//-------------------------------------------------------------------// 
//  Synopsis:   Information about poles of H(s)=num(s)/den(s)
//  Syntax:     damp ( num ) 
//  Description: 
//  
//  Example:
//  H=(num,den);
//  num=[5,3,1];
//  den=[1,6,4,4];
//  D=damp(den);
//  D.zeta=
//  D.wn=
//  D.tau=
//  D.Q=
//  D.P=
//  D.Fpeak=
//  D.overshoot=
//  
//  Example 2:
// RLC
// L=15.9e-3;
// C=.159e-3;
// R=10;
// R=13;
// R=13;
// num2=[1/L,0];
// den2=[1,R/L,1/(L*C)];

// X=damp(den2);

//        for(k in members(X)){printf("%s\r\t\t = %f\n",k,X.[k]);}
//           FRpeaking        = 0.105613
//           Fringing         = 76.067475
//           P                = -408.805031
//           Q                = 0.384615
//           overshoot        = 0.068077
//           tau              = 0.002446
//           wn               = 628.930818
//           zeta             = 0.650000


//-------------------------------------------------------------------//
require angle roots freqresp_peaking overshoot

damp  = function(den, print)
{
global(pi)

# Damping Ratio
zeta=-cos(angle(roots(den)));
# Another formula ofDamping Ratio
# zeta=-real(roots(den))./wn;

# Natural frequency
wn=abs(roots(den));

# Time constant
tau=2/(wn.*zeta);

# Pole location
P=roots(den);

Q=1/(2*zeta);

if(max(size(den)==3)){FRpeaking=freqresp_peaking(zeta); }
if(max(size(den)==3)){OS=overshoot(zeta); }

# Ringing frequency
if(max(size(den)==3)){ Fringing = wn./(2*pi) .* sqrt(1- zeta[1]*zeta[1]);}


# Logarithic decrement is defined as the natural log of the ratio of any two successive amplitudes
# delta = ln(x1/x2)
# zeta = delta/(2*pi)

X =  << zeta=zeta; wn=wn; tau=tau; P=P; Q=Q; overshoot=OS; FRpeaking=FRpeaking; Fringing=Fringing  >>;

# if(nargs == 2){
	# for(k in members(X)){printf("%s\r\t\t = %f\n",k,X.[k]);}
# 
# }

return X;
  // return atan2(imag(a), real(a));
};
