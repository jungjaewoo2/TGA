# EXAMPLE 1:

require tfd
require rem

#
# Make a chirp
#

x = exp(2*pi*1j* 6*((1:40)/40).^2);
show(x);
plot (x.');
pause("Hit return for contour");

fa=20; sigma=inf();
</ max_t; min_t; y /> = tfd (x,fa,sigma);

plcont (<< x=1:y.nr; y=1:y.nc; z=y >>);
pause("Hit return for instantaneous frequency");

#
# Plot instantaneous frequency
#

plegend();
pltitle ("Instanteous Frequency");
plot( [(min_t:max_t)', (12*(min_t:max_t)/40)' ] );
pause("Hit return for contour");

#
# Plot mesh
#

pltitle ("Wigner-Ville Transform");
plmesh (<< x=1:y.nr; y=1:y.nc; z=y >>);
