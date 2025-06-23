//function [iu,nargs,y,y2]=mulresp(fun,a,b,c,d,t,nargo,bodeflag);
//MULRESP Multivariable response.
//
//	[IU,NARGS,Y,Y2] = MULRESP('fun',A,B,C,D,T,NARGO,BODEFLAG)


mulresp=function(fun,a,b,c,d,t,nargo,bodeflag)
{

  r=d.nr;
  m=d.nc;
  if (r*m>1) {
    // MIMO system
    iu=0;
    if (nargo==0) {
	clg();
	hold off
	if (r*m==2) { sp=210; else sp=220; }
	if (r*m>4||bodeflag==1) {
	   disp("Strike any key after each screen");
	}
	scnt=0;
	for (i in 1:m) {
		if (bodeflag==0) {
		  for (j in 1:r) {
		    if (scnt==4) {
		        pause(); clg(); scnt=0;
		    }
		    scnt=scnt+1; 
		    subplot(sp+scnt)
                    if (!isempty(c)) { cj = c[j;]; }
		    if (!isempty(d)) { dj = d[j;]; }
		    fun(a,b,cj,dj,i,t);
		    pltitle("Input "+int2str(i)+" Output "+int2str(j));
		  }
		else 
		  if (scnt==4) { pause(); clg(); scnt=0; }
		  scnt = scnt+4;
		  fun(a,b,c,d,i,t);
		  pltitle("Input "+int2str(i));
		}
	}
	subplot(111)
    else
	y=[]; 
	y2=[];
	for (i in 1:m) {
		if (bodeflag==0) {
		  y=[y,fun(a,b,c,d,i,t)];
		else
		  // Force compile to recognize these variables: 
		  phase =[];  mag = [];
		  eval(["[mag,phase]=",fun,"(a,b,c,d,i,t);"])
		  y  = [y,mag]; 
		  y2 = [y2,phase];
		}
	}
    }
  else		
    // SISO systems
    iu=1; 
    nargs=5;
  }

  return <<iu=iu; nargs=nargs; y=y; y2=y2>>;	
};

