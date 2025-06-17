dfdp = function(x,p){
	global(pi)
	// jacobian df/dp
      y=[];
     y[;1] = (p[3].*x.*exp(-x/p[1]).*sin(2*pi*p[2]*x))/(p[1].*p[1]);
     y[;2] = (2*pi*p[3].*x.*exp(-x/p[1]) .* cos(2*pi*p[2].*x));
     y[;3] = (exp(-x/p[1]).* sin(2*pi*p[2].*x));
     y[;4] = 1;
     y[;5] =-2*p[2]* p[3]*pi*exp((-x)./p[1]) .* cos(2*p[2]*pi*(x-p[5]));

return y;
};

