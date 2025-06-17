tfeval = function(ns,ds,f){
 //   This function evaluates a transfer function for a given frequency(s)
 //
 //       p = tfeval(ns, ds, f)
 //
 //   where NS => the vector representing the numerator polynomial
 //        DS => the vector representing the denominator polynomial
 //        F  => vector of frequency values in Hz at which the TF is
 //              to be evaluated
 //        P  => vector of TF values corresponding to F.
 //
 //                ns(1)*s^N + ns(2)*s^(N-1) ...+ ns(N)*s^1 + ns(N+1)
 //       H(s) = ---------------------------------------------------
 //                ds(1)*s^M + ds(2)*s^(M-1) ...+ ds(M)*s^1 + ds(M+1)
 //
 //     Example:
 //               for linear frequency evaluation
 //               >> freq = [0:10:10000];   %  Frequency points
 //               >> numc = [10  0];   %  Numerator coefficients
 //               >> denc = [1 5000*2*pi];   %  Denominator coefficients 
 //               >> tf = tfeval(numc, denc, freq);    % Compute complex TF points
 //               >>  figure; plot(f,abs(p)); ylabel('Magnitude')   %  Plot magnitude
 //               >>  figure; plot(f,angle(p)*180/pi); ylabel('Degrees')  % Plot Phase
 //
 //       Updated by Kevin D. Donohue  ( donohue@engr.uky.edu ) February 1, 2010
 
 global(pi)
 num_order=length(ns)-1; //  Determine order of numerator
 den_order=length(ds)-1; //  Determine order of denominator
 
 s=2i*pi*f; //  Create s vector for evaluating TF
 
 //  Loop to sum up every term in numerator 
 sumn=zeros((s)); // Initialize accumulation variable
 for (k in 1:num_order+1){  
   sumn=sumn+ns[k]*s.^(num_order-k+1); 
   }
 
 //  Loop to sum up every term in numerator 
 sumd=zeros((s)); // Initialize accumulation variable
 for (k in 1:den_order+1){  
   sumd=sumd+ds[k]*s.^(den_order-k+1); 
   }
 p=sumn./sumd; //  Divide numerator and denominator points
 
 
 
 return p; 
 };
