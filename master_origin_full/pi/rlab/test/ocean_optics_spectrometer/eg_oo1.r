//
//
//
gnuwins(1);
clear (savefile);

rfile libocean.so

// initialize library and assign spectrometer No.0 to variable 'spec'
spec = ocean_init(0);

// set integration time
spec.int_time(1e-2);

// no boxing average of signal
spec.boxcar_width(0);

// average two scans
spec.scan_avg(2);

// correct for dark current
spec.corr_dark(1);

data =<<>>;
for (i in 1:100)
{
	spinner();
	ix = text(i,"%03.0f");

	data.[ix] = spec.get_spectrum();

	rfile module_plot

	sleep(30);
}

savefile="./spectrum.eps";
rfile module_plot
