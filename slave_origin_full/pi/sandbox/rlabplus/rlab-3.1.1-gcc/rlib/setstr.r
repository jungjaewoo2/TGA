//-------------------------------------------------------------------//

// Synopsis:	Translates a vector of integers to a string of
//		ASCII characters.

// Syntax:	setstr ( n )

// Description:

//	Setsrt converts a vector of integers in n to a string of ASCII
//	characters.  Only integer values between 0 and 127 will be
//	converted successfully.  Nonprintable characters will be
//	translated to their symbols surrounded by braces.

//	This is useful when fread is used to read "char" type data
//	from a file.  Fread returns an integer vector.

//	name=fread(file_name,8,"char");
//	sname=setstr(name);

//	The example lines above will read in 8 bytes of characters into
//	integer vector name.  Then name is converted to a string sname
//	by setstr.

//	This file is a substitute for the same in matlab by Derrick Early.

//-------------------------------------------------------------------//

setstr = function(n)
{
  _ascii=["{nul}","{soh}","{stx}","{etx}","{eot}","{enq}","{ack}","{bel}", ...
         "{bs}", "{ht}", "{nl}", "{vt}", "{np}", "{cr}", "{so}", "{si}",  ...
         "{dle}","{dc1}","{dc2}","{dc3}","{dc4}","{nak}","{syn}","{etb}", ...
         "{can}","{em}", "{sub}","{esc}","{fs}", "{gs}", "{rs}", "{us}",  ...
         " ",    "!",    "\"",   "#",    "$",    "%",    "&",    "'",     ...
         "(",    ")",    "*",    "+",    "\,",   "-",    ".",    "/",     ...
         "0",    "1",    "2",    "3",    "4",    "5",    "6",    "7",     ...
         "8",    "9",    ":",    ";",    "<",    "=",    ">",    "?",     ...
         "@",    "A",    "B",    "C",    "D",    "E",    "F",    "G",     ...
         "H",    "I",    "J",    "K",    "L",    "M",    "N",    "O",     ...
         "P",    "Q",    "R",    "S",    "T",    "U",    "V",    "W",     ...
         "X",    "Y",    "Z",    "[",    "\\",   "]",    "^",    "_",     ...
         "`",    "a",    "b",    "c",    "d",    "e",    "f",    "g",     ...
         "h",    "i",    "j",    "k",    "l",    "m",    "n",    "o",     ...
         "p",    "q",    "r",    "s",    "t",    "u",    "v",    "w",     ...
         "x",    "y",    "z",    "{",    "|",    "}",    "~",    "{del}"];
  S="";

  for (i in 1:length(n))
  { S=S+_ascii[n[i]+1]; }

  return S;
};
