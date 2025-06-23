//
// file: eg_writer1.r
//
// based on writer-example1 from perl package ooolib:
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
// User defined variables
ooo.set ("meta1-name", "Programmers");                         // 1-4 user defined variables
ooo.set ("meta1-value", "Joseph Colton and Marijan Kostrun");  // 1-4 user defined variables


for (i in 1:9)
{ ooo.data("h"+text(i), "Heading "+text(i)); }

ooo.data("textbody", ...
    [ "The lazy fox jumped over the laying log. The log sat under the jumping fox.\n"+ ...
      "The lazy fox jumped over the laying log. The log sat under the jumping fox.\n"+ ...
      "The lazy fox jumped over the laying log. The log sat under the jumping fox."] );

ooo.generate("writer-example1.sxw");

