//
//
//
gnuwins(1);
clear (savefile);

rfile libocean.so

// initialize library and assign spectrometer No.0 to variable 'spec'
spec = ocean_init(0);

// set as integration time the minimum time allowed by spectrometer
spec.int_time(1e-2);

// no boxing of signal
spec.boxcar_width(0);

// average two scans
spec.scan_avg(2);

// correct for dark current
spec.corr_dark(1);

data = <<>>;
tic();
data.test = spec.get_spectrum();
toc()

rfile module_plot
