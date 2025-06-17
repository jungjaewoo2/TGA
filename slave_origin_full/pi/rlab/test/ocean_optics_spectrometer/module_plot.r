//
//
//

gnuwin (1);
gnulimits(400,800);
gnuxlabel("Wavelength (nm)");
gnuxtics (50,5);
gnuylabel("Intensity");
gnulegend(members(data));
if (exist(savefile))
{
	gnuplot (data, savefile);
else
	gnuplot (data);
}

