// based on writer-example5 from perl package ooolib:
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

for (s in [6:96:2])
{
  r = s * 2;
  b = 256 - r;
  c = text(int(r), "%2x") + "99" + text(int(b),"%2x");
  ooo.set ("text-color", c);
  ooo.set ("text-size",  s);
  ooo.data("default", "ABC");
}

ooo.generate("writer-example5.sxw");

