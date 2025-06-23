//
// file: eg_calc2.r
//
// based on calc-example2 from perl package ooolib:
//  This program creates a SXC Calc file with 25 cells and a formula.
//  by Joseph Colton

rfile libooo

# Set variables
ooo.init("sxc");
ooo.set ("builddir", "./tmp");
ooo.set ("title", "Calc Column Sizes");
ooo.set ("subject", "Calc Example");
ooo.set ("comments", "Sizes are set using cm");
ooo.set ("author", "ooolib Example");
// User defined variables
ooo.set ("meta1-name", "Programmers");                         // 1-4 user defined variables
ooo.set ("meta1-value", "Joseph Colton and Marijan Kostrun");  // 1-4 user defined variables

// Write the spreadsheet data
// It is in x, y cords.
for (x in 1:5)
{
  ooo.set ("column-width", x, x);
  for(y in 1:5)
  {
    ooo.set ("cell-loc", x, y);
    ooo.data("cell-float", x*y);
  }
}
ooo.set ("cell-loc", 1 ,6);
ooo.data("cell-formula", "=SUM(B2:C2)");
ooo.generate("calc-example2_rlab.sxc"); // Create the document

