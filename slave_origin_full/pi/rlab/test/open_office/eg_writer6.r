// based on writer-example6 from perl package ooolib:
//  This program creates a SXW writer file
//  by Joseph Colton

rfile libooo

// Set variables
ooo.init("sxw");
ooo.set ("builddir", "./tmp");                   // Directory to create document in
ooo.set ("title", "OpenOffice.org Writer Document");
ooo.set ("subject", "Writer Examples");
ooo.set ("comments", "Options can be turned on, then off again.");
ooo.set ("keyword", "keyword1");                 // Add a keyword for the document
ooo.set ("keyword", "keyword2");                 // Add a keyword for the document
ooo.set ("author", "libooo example");            // Set the author of the document

ooo.data("h1", "Tropical Fruit");

lst = ooo.special("list-ordered", "new", "Pineapple");
ooo.special("list-ordered", lst, "Papaya");
ooo.special("list-ordered", lst, "Banana");
ooo.special("list-ordered", lst, "Mango");

ooo.set ("text-font", "Bookman");
ooo.data("textbody", "Bookman");
ooo.set ("text-font", "Courier 10 Pitch");
ooo.data("textbody", "Courier 10 Pitch");
ooo.set ("text-font", "Zapf Chancery");
ooo.data("textbody", "Zapf Chancery");

ooo.generate("writer-example6.sxw");

