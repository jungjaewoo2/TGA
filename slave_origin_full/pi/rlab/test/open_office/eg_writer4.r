// based on writer-example4 from perl package ooolib:
//  This program creates a SXW writer file
//  by Joseph Colton

rfile libooo

// Set variables
ooo.init("sxw");
ooo.set ("builddir", "./tmp");                   // Directory to create document in
ooo.set("title", "OpenOffice.org Writer Document");
ooo.set("subject", "Writer Examples");
ooo.set("comments", "Options can be turned on, then off again.");
ooo.set ("keyword", "keyword1");                 // Add a keyword for the document
ooo.set ("keyword", "keyword2");                 // Add a keyword for the document
ooo.set ("author", "libooo example");            // Set the author of the document

ooo.set("text-color", "blue");
ooo.data("default", "blue text");
ooo.set("text-color", "red");
ooo.data("default", "red text");
ooo.set("text-color", "green");
ooo.data("default", "green text");

ooo.set("text-color", "default");
ooo.data("default", "Be prepared for visual pain.");
ooo.special("pagebreak");

ooo.set("text-bgcolor", "red");
ooo.set("text-color", "blue");
ooo.data("default", "blue on red text");
ooo.set("text-bgcolor", "green");
ooo.set("text-color", "red");
ooo.data("default", "red on green text");
ooo.set("text-bgcolor", "blue");
ooo.set("text-color", "green");
ooo.data("default", "green on blue text");

ooo.generate("writer-example4.sxw");

