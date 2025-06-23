//
// file: eg_calc4.r
//
// based on calc-example4 from perl package ooolib:
//  This program creates a SXC Calc file
//  by Joseph Colton

rfile libooo

// Set variables
ooo.init("sxc");
ooo.set ("builddir", "./tmp");                   // Directory to create document in
ooo.set ("title", "Calc 25 Cells");              // Title of document
ooo.set ("subject", "Calc Example");             // Subject of document
ooo.set ("comments", "Who needs comments?");     // Comments for document
ooo.set ("keyword", "keyword1");                 // Add a keyword for the document
ooo.set ("keyword", "keyword2");                 // Add a keyword for the document
ooo.set ("author", "libooo example");            // Set the author of the document

// User defined variables
ooo.set ("meta1-name", "Programmers");                         // 1-4 user defined variables
ooo.set ("meta1-value", "Joseph Colton and Marijan Kostrun");  // 1-4 user defined variables

// Write the spreadsheet data
ooo.data("cell-float", 1);
ooo.set ("cell-down");
ooo.data("cell-float", 2);
ooo.set ("cell-down");
ooo.data("cell-float", 3);
ooo.set ("cell-right");
ooo.data("cell-float", 4);
ooo.set ("cell-right");
ooo.data("cell-float", 5);
ooo.set ("cell-up");
ooo.data("cell-float", 6);
ooo.set ("cell-up");
ooo.data("cell-float", 7);
ooo.set ("cell-left");
ooo.data("cell-float", 8);

ooo.generate("calc-example4.sxc");

