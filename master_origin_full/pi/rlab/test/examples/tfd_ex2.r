# EXAMPLE 2:

require tfd
require rem

#
# Two crossing chirps...
#

x = exp(2*pi*1j*6*((1:80)/80).^2) + ...
    exp(2*pi*1j*(1:80)/80.*(50-30*(1:80)/80));
plegend();
plot (x.');
pause ("Hit return to see WVD distribution");

fa = 40;
sigma = inf();

#
# Make a WVD
#

</ max_t ; min_t ; y /> = tfd (x,fa,sigma);

plcont (<< x=1:y.nr; y=1:y.nc; z=y >>);
pause("Hit return to see WVD");

plmesh (<< x=1:y.nr; y=1:y.nc; z=y >>);
pause("Hit return to see CWD");

#
# Make a CVD
#

sigma = 1;
</ max_t ; min_t ; y /> = tfd(x,fa,sigma);

plcont (<< x=1:y.nr; y=1:y.nc; z=y >>);

pause("Hit return to see CVD Transform");
plmesh (<< x=1:y.nr; y=1:y.nc; z=y >>);

pause("Hit return to see Instaneous frequency");

plot (<< [(min_t:max_t)',(12*(min_t:max_t)/80)'] ;
         [(min_t:max_t)',(50-60*(min_t:max_t)/80)'] >>);
