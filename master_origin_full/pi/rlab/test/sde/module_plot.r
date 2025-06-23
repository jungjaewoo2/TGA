//
// plot brownian path
//

havegnuwin(1,,"test.gnu");

gnucmd   (1,"set title 'Strong convergence of E-M method'");
gnuxlabel("Time step");
gnuylabel("Error");
gnulimits(0.001,0.1, 0.1,1);
gnuytics ([0.1,0.2,0.5,1]);
gnuxtics ([0.001,0.002,0.005,0.01, 0.02, 0.05, 0.1]);
gnucmd   (1, "unset grid; set grid xtics ytics;");
gnuscale ("log", "log");
gnuformat(["with lines"]);
gnuplot  (1, data);
