//----------------------------------------------------------------------
// parallel
//
// syntax:  </a;b;c;d/> = parallel(a1,b1,c1,d1,a2,b2,c2,d2,e,f,g,h)
//          </den;num/> = parallel(num1,den1,num2,den2) 
//
//	Parallel connection of two systems.  
//
//                +-->[System1]--+
//            u-->+              O--->y
//                +-->[System2]--+
//
//	</A;B;C;D/> = parallel(A1,B1,C1,D1,A2,B2,C2,D2)  produces a state-
//	space system consisting of the parallel connection of systems 1 
//	and 2 that connects all the inputs together and sums all the 
//	outputs of the two systems,  Y = Y1 + Y2.
//
//	</A;B;C;D/> = parallel(A1,B1,C1,D1,A2,B2,C2,D2,IN1,IN2,OUT1,OUT2) 
//	connects the two systems in parallel such that the inputs 
//	specified by IN1 and IN2 are connected and the outputs specified
//	by OUT1 and OUT2 are summed. The vectors IN1 and IN2 contain 
//	indexes into the input vectors of system 1 and system 2, 
//	respectively.  Similarly, the vectors OUT1 and OUT2 contain 
//	indexes into the outputs of the systems.  The parallel connection
//	is performed by appending the two systems, summing the specified
//	inputs and outputs, and removing the, now redundant, inputs and 
//	outputs of system 2.
//
//	</DEN;NUM/> = parallel(NUM1,DEN1,NUM2,DEN2) produces a parallel 
//	connection of the two transfer function systems.
//
//	See also: cloop, feedback, and series. 
//
//----------------------------------------------------------------------
require abcdchk append ssdelete tfchk

parallel=function(a1,b1,c1,d1,a2,b2,c2,d2,e,f,g,h)
{

  // --- Determine which syntax is being used ---
  if (nargs == 4) {  
    // Form Parallel connection of T.F. system ---
    tmp1 = tfchk(a1,b1); 
    tmp2 = tfchk(c1,d1);
    nn = tmp1.numc.nr;
    mn = tmp1.numc.nc;
    for (k in 1:nn) {
      a[k;] = conv(tmp1.numc[k;],tmp2.denc) + conv(tmp2.numc[k;],tmp1.denc);
      b = conv(tmp1.denc,tmp2.denc);
    }
    return <<num=a;den=b>>;
    
  } else { if (((nargs>=5) && (nargs<=7)) || ((nargs>=9) && (nargs<=11))) {
    error("Wrong number of input arguments.");

  } else { if ((nargs==8) || (nargs==12)) {
    // State space systems 
    msg = abcdchk(a1,b1,c1,d1);
    if (msg != "") { error("system 1:"+msg); }
    msg = abcdchk(a2,b2,c2,d2);
    if (msg != "") { error("system 2:"+msg); }
    ny1 = d1.nr; nu1 = d1.nc;
    ny2 = d2.nr; nu2 = d2.nc;
    if (nargs == 8) {
      // State space systems w/o selection vectors
      inputs1 = [1:nu1];     outputs1 = [1:ny1];
      inputs2 = [1:nu2]+nu1; outputs2 = [1:ny2]+ny1; 
    }
    if (nargs == 12) {
      // State space systems with selection vectors
      inputs1 = e;      outputs1 = g;
      inputs2 = f+nu1;  outputs2 = h+ny1;
    }
  
    // Check sizes
    if (length(inputs1)!=length(inputs2)) {
       error("Input sizes don't match.");
    }
    if (length(outputs1)!=length(outputs2)) {
       error("Output size don't match."); 
    }

    // --- Parallel Connection ---
    </a;b;c;d/> = append(a1,b1,c1,d1,a2,b2,c2,d2);

    // Connect inputs
    if (!isempty(b)) { b[;inputs1]=b[;inputs1]+b[;inputs2]; }
    if (!isempty(d)) { d[;inputs1]=d[;inputs1]+d[;inputs2]; }
 
    // Connect outputs
    if (!isempty(c)) { c[outputs1;]=c[outputs1;]+c[outputs2;]; }
    if (!isempty(d)) { d[outputs1;]=d[outputs1;]+d[outputs2;]; }

    </a;b;c;d/> = ssdelete(a,b,c,d,inputs2,outputs2);
    return <<a=a;b=b;c=c;d=d>>;
  
  } else {
    error("Wrong number of input argumets.");
    
  }}}

};

