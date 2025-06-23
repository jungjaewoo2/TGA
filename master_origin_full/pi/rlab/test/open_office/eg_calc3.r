//
// file: eg_calc1.r
//
// based on calc-example1 from perl package ooolib:
//  This program creates a SXC Calc file with 25 cells and a formula.
//  by Joseph Colton

rfile libooo

// Set variables
ooo.init("sxc");
ooo.set ("builddir", "./tmp");                       // Directory to create document in
ooo.set ("title", "Calc 25 Cells");              // Title of document
ooo.set ("subject", "Calc Example");             // Subject of document
ooo.set ("comments", "Who needs comments?");     // Comments for document
ooo.set ("keyword", "keyword1");                 // Add a keyword for the document
ooo.set ("keyword", "keyword2");                 // Add a keyword for the document
ooo.set ("author", "libooo example");            // Set the author of the document
ooo.set ("cell-auto", 0, 1);

// User defined variables
ooo.set ("meta1-name", "Programmers");                         // 1-4 user defined variables
ooo.set ("meta1-value", "Joseph Colton and Marijan Kostrun");  // 1-4 user defined variables

# Write the spreadsheet data
for (x in 1:10)
{ ooo.data("cell-float", x); }

ooo.set ("cell-loc", "C3");
ooo.data("cell-float", "123.697699");

ooo.generate("calc-example3_rlab.sxc"); // Create the document

