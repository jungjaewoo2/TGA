ResRise = function (Tau_decay, fres)
{
// Convert resonant rise to decay
// Convert decay to resonant rise
//
// STATUS: NEW
//
// Arguments:
//      Tau_decay          - decay time in [ms]
//      fres               - resonant frequency in [Hz]
// Return:
//      Resonant Rise      - in [dB]
//
//
// Example:

//            fres = [93.5698; ...
//                    92.2385; ...
//                    92.9085; ...
//                    96.9862; ...
//                    96.4145];
 
 
//            decay = [3.96085; ...
//                     4.05202; ...
//                     4.12316; ...
//                     3.76456; ...
//                     3.69673];


//   ResRise (decay, fres) = 
//                           1.32147939  
//                           1.39467369  
//                           1.60871018  
//                           1.19148079  
//                           0.98219916 


global(pi)

return 20* log10(pi .* Tau_decay .* 1e-3.* fres)
};  
