verify_wiring = function (fn) {
// fn -> filename or matrix with fine motion data (includes laser analog output)
//
// verifies polarity of winding (and correct magnet installation
//
//
// Usage:  verify_wiring("/tmp/prod.dat")
//         or
//         verify_wiring(d)
//
// Output: +1 normal polarity
//         -1 reversed polarity (miswired coil)
//
/////////////////////////////////////////////////////////////////////////////////

global(d)

if(class(fn)== "string"){if(isfile(fn)){ xx=read_ascii(fn);}}

if(class(fn)== "num"){d=fn;}



d0=d-(mean(d));
pol= sign(sum(d0[;3].*d0[;5]));


return pol;
};
