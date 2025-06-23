// based on writer-example3 from perl package ooolib:
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

// Left side
ooo.data("default", "Regular left");
ooo.set ("bold", "on");
ooo.data("default", "Bold left");
ooo.set ("bold", "off");
ooo.set ("italic", "on");
ooo.data("default", "Italic left");
ooo.set ("italic", "off");
ooo.set ("underline", "on");
ooo.data("default", "Underline left");
ooo.set ("underline", "off");

// Center
ooo.set ("justify", "center");
ooo.data("default", "Regular center");
ooo.set ("bold", "on");
ooo.data("default", "Bold center");
ooo.set ("bold", "off");
ooo.set ("italic", "on");
ooo.data("default", "Italic center");
ooo.set ("italic", "off");
ooo.set ("underline", "on");
ooo.data("default", "Underline center");
ooo.set ("underline", "off");

// Right
ooo.set("justify", "right");
ooo.data("default", "Regular right");
ooo.set("bold", "on");
ooo.data("default", "Bold right");
ooo.set("bold", "off");
ooo.set("italic", "on");
ooo.data("default", "Italic right");
ooo.set("italic", "off");
ooo.set("underline", "on");
ooo.data("default", "Underline right");
ooo.set("underline", "off");

// Block
ooo.set("justify", "block");
ooo.data("default", "Regular block");
ooo.set("bold", "on");
ooo.data("default", "Bold block");
ooo.set("bold", "off");
ooo.set("italic", "on");
ooo.data("default", "Italic block");
ooo.set("italic", "off");
ooo.set("underline", "on");
ooo.data("default", "Underline block");
ooo.set("underline", "off");

ooo.generate("writer-example3.sxw");

