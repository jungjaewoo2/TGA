//
// file: eg_calc5.r
//
// based on calc-example5 from perl package ooolib:
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

// Column 1
ooo.set ("bold", "on");
ooo.set ("cell-loc", 1, 1);
ooo.data("cell-text", "bold");

ooo.set ("text-color", "ff0000");
ooo.set ("cell-loc", 1, 2);
ooo.data("cell-text", "red");

ooo.set ("text-color", "00ff00");
ooo.set ("cell-loc", 1, 3);
ooo.data("cell-text", "green");

ooo.set ("text-color", "0000ff");
ooo.set ("cell-loc", 1, 4);
ooo.data("cell-text", "blue");
ooo.set ("bold", "off");
ooo.set ("text-color", "000000");

// Column 2
ooo.set ("italic", "on");
ooo.set ("cell-loc", 2, 1);
ooo.data("cell-text", "italic");

ooo.set ("text-color", "ff0000");
ooo.set ("cell-loc", 2, 2);
ooo.data("cell-text", "red");

ooo.set ("text-color", "00ff00");
ooo.set ("cell-loc", 2, 3);
ooo.data("cell-text", "green");

ooo.set ("text-color", "0000ff");
ooo.set ("cell-loc", 2, 4);
ooo.data("cell-text", "blue");
ooo.set ("italic", "off");
ooo.set ("text-color", "000000");

ooo.generate("calc-example5.sxc");
