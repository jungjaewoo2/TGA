// based on writer-example2 from perl package ooolib:
//  This program creates a SXW writer file
//  by Joseph Colton

rfile libooo

// Set variables
ooo.init("sxw");
ooo.set ("builddir", "./tmp");                   // Directory to create document in
ooo.set("title", "OpenOffice.org Writer Document");
ooo.set("subject", "Writer Examples");
ooo.set("comments", "Simple Writer Document");
ooo.set ("keyword", "keyword1");                 // Add a keyword for the document
ooo.set ("keyword", "keyword2");                 // Add a keyword for the document
ooo.set ("author", "libooo example");            // Set the author of the document

ooo.data("default", "This is some simple text");
ooo.data("default"); // A blank line
ooo.data("default", "Special Characters: ");
ooo.data("default", "\t!\"#\$%\&'()`{}[];:+*\\_,.<>?|");

ooo.generate("writer-example2.sxw");

