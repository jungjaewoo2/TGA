//----------------------------------------------------------------------
//
// printsys
//
// syntax: printsys(A,B,C,D)
//         printsys(A,B,C,D,ULABELS,YLABELS,XLABELS)
//         printsys(NUM,DEN,"s") 
//         printsys(NUM,DEN,"z")
//
//      Print system in pretty format.
//
//	printsys is used to print state space systems with labels to the
//	right and above the system matrices or to print transfer functions
//	as a ratio of two polynomials.
//
//	printsys(A,B,C,D,ULABELS,YLABELS,XLABELS) prints the state-space
//	system with the input, output and state labels contained in the
//	strings ULABELS, YLABELS, and XLABELS, respectively.  ULABELS, 
//	YLABELS, and XLABELS are string variables that contain the input,
//	output and state labels delimited by spaces.  For example the 
//	label string
//		YLABELS=["Phi", "Theta", "Psi"]
//	defines "Phi" as the label for the first output, "Theta" as the 
//	label for the second output, and "Psi" as the label for the third
//	output.
//	
//	printsys(A,B,C,D) prints the system with numerical labels.
//	
//	printsys(NUM,DEN,"s") or printsys(NUM,DEN,"z") prints the transfer
//	function as a ratio of two polynomials in the transform variable
//	"s" or "z".
//
//	See also: printmat and poly2str.
//
//----------------------------------------------------------------------
require abcdchk isstr poly2str printmat tfchk

printsys = function(a,b,c,d,ulab,ylab,xlab)
{

  len = 12; // Label parameter.  Should match parameter in PRINTMAT.

  // --- Determine which syntax is being used ---
  if (nargs==2) {		
    // Transfer function without transform variable
    T = tfchk(a,b);
    num = T.numc;
    den = T.denc;
    tvar = "s";
  else if (nargs==3) {	
    // Transfer function with transform variable
    T = tfchk(a,b);
    num = T.numc;
    den = T.denc;
    if (!isstr(c)) { error("Transform variable must be a string."); }
    tvar = c;
  else if (nargs==4) {
    // State space system without labels
    msg = abcdchk(a,b,c,d);
    if (msg!="") { error("printsys:"+msg); }
    nx = a.nr;
    na = a.nc;
    if (nx>0) { 
      ny = c.nr;
      nc = c.nc;
      nb = b.nr;
      nu = b.nc;
    else
      ny = d.nr;
      nu = d.nc;
    }
    ulab = []; 
    ylab = []; 
    xlab = [];
    for (i in 1:nx) {
      sprintf(T,"        x%-3.0f",i);
      xlab = [xlab;T];
    }
    printmat(a,"A",xlab,xlab);
    for (i in 1:nu) {
      sprintf(T,"        u%-3.0f",i);
      ulab = [ulab; T];
    }
    printmat(b,"B",xlab,ulab);
    for (i in 1:ny) {
      sprintf(T,"        y%-3.0f",i);
      ylab = [ylab; T];
    }
    printmat(c,"C",ylab,xlab);
    printmat(d,"D",ylab,ulab);
  else if (nargs==5||nargs==6) {
    error("Wrong number of input arguments.");
  else			
    // State space system with labels
    msg = abcdchk(a,b,c,d);
    if (msg!="") { error("printsys:"+msg); }

    nx = a.nr;
    na = a.nc;
    ny = c.nr;
    nc = c.nc;
    nb = b.nr;
    nu = b.nc;
    if (!isstr(ulab)) { error("ULAB must be a string."); }
    if (!isstr(ylab)) { error("YLAB must be a string."); }
    if (!isstr(xlab)) { error("XLAB must be a string."); }
  
    // Put labels into row-wise format for MATPRINT
    printmat(a,"A",xlab,xlab);
    printmat(b,"B",xlab,ulab);
    printmat(c,"C",ylab,xlab);
    printmat(d,"D",ylab,ulab);
  }}}}

  if (nargs==2||nargs==3){
    // Print system as a ratio of polynomials
    nn = num.nr;
    mn = num.nc;
    T = poly2str(den,tvar);
    s = T.sout;
    dlen = T.lenout;
    for (i in 1:nn) {
      T = poly2str(num[i;],tvar);
      t = T.sout;
      nlen = T.lenout;
      // Now print the polynomials
      len = max([dlen,nlen])-3;
    
      if (nn==1) {
        printf("\nnum/den = \n");
      else
        printf("\nnum[%i]/den = \n",i);
      }
    
      if (length(t)<len+3) {
        // Print numerator
        for (i in 1:(len+4-length(t))/2) { printf(" "); } printf("%s\n",t);
      else
        printf("%s\n",t);
      }
      printf("   "); for (i in 1:len) { printf("-"); } printf("\n");
      if (length(s)<len+3) {
        // Print denominator
        for (i in 1:(len+4-length(s))/2) { printf(" "); } printf("%s\n",s);
      else
        printf("%s\n",s);
      }
    }

  }
};
