//   writeltx.r
//   J.J.Green
//   Version 0.9 06/03/97
//   This code is in the public domain.

//---------------------------------------------------------------------------//

// Synopsis:	Write a matrix in LaTeX format.

// Syntax:	writeltx(<filename>,<matrix>,<options>)

// Description:

//   Writes <matrix> to the file <filename> in the LaTeX format.
//
//   The optional third argument is a list with the following members:
//
//     dp       : decimal places                   (5)
//     width    : width                            (8)
//     format   : exponential, float or integer    ("float")
//     imUnit   : LaTeX code to be used for the
//                imaginary unit                   ("i")
//     align    : alignment of matrix columns      ("centre")
//
//   It is important to note that Rlab uses C-style escape sequences
//   (such as \n = newline, \t = tab) in strings. So the LaTeX backslash
//   must be written as \\ in the imUnit string. For example if one wanted
//   a bold-face j as the imaginary unit then something like
//
//     opt.imUnit="\\textbf{j}";
//     writeltx("file.tex",a,opts);
//
//   would be required.
//
//   Note also that the first call to writeltx() overwrites <filename>
//   and subsequent calls are appended. You need to close the file
//   formally (as with writem) when you have finished using it.
//
//   To do: Testing, support scientific exponential notation (x10^3n)
//
//   Reference:
//     H.Kopka & P.W.Daly,`A Guide to LaTeX2e' (second edition),
//     Addison-Wesley, 1995.
//

// Dependencies: none

// the following are internal functions so declared static
static(_formatReal,_formatComplex)

writeltx = function(filename, matrix, options){

  // Declare local variables
  local(rows,columns,columnSeparator,rowLoop,columnLoop);
  local(rowString,numberString,alignchar);

  // Check that the arguments are feasible
  if (!((nargs == 2)||(nargs == 3))){
    printf("Usage: writeltx(<filename>,<matrix>,<options>)\n");
    printf("see help(write) for list of options\n");
    return(0);
    }
  if (! exist(filename)){
    printf("First argument does not exist\n");
    return(0);
    }
  if (! exist(matrix)){
    printf("Second argument does not exist\n");
    return(0);
    }
  if ( type(filename) != "string"){
    printf("First argument must be a string\n");
    return(0);
    }
  if ( !((type(matrix) == "real") || (type(matrix) == "complex"))){
    printf("Second argument must be a real or complex matrix\n");
    return(0);
    }

  // Set default options
  if (! exist(options)){options=<<>>;}
  if (! exist(options.dp)){options.dp=5;}
  if (! exist(options.width)){options.width=8;}
  if (! exist(options.format)){options.format="float";}
  if (! exist(options.align)){options.align="centre";}
  if (! exist(options.imUnit)){options.imUnit="i";}

  // make the format string
  if (options.format=="float"){
    sprintf(options.fstring,"%%%i.%if",options.width,options.dp);
    } else { {
      if (options.format=="exponential"){
        sprintf(options.fstring,"%%%i.%ie",options.width,options.dp);
        } else { {
          if (options.format=="integer"){
            sprintf(options.fstring,"%%%i.%ii",options.width,options.dp);
            } else { {
              printf("Unknown format type \"%s\"\n",options.format);
              return(0);
              }
            }
          }
        }
      }
    }

  // How do we align columns in the matrix? -- just get the first
  // character of the options.align string (thus we avoid trans-Atlantic
  // arguments about `center' versus `centre')
  alignchar=strsplt(options.align)[1];

  // Set some constants
  rows = matrix.nr;
  columns = matrix.nc;
  columnSeperator = " & ";

  // Open the %%\LaTeX\ %% file,
  if (close(filename)){
    // if it was open then we append...
    open(filename,"a");
    } else {{
      // ... otherwise overwrite it.
      open(filename,"w");
      }
    }

  // Print the preamble
  fprintf(filename,"%% RLab writeltx() output\n");
  fprintf(filename,"\\begin{equation}%%\\label{}\n");
  fprintf(filename,"\\left[\\begin{array}{");
  for (columnLoop in 1:columns){
    fprintf(filename,"%s",alignchar);
    }
  fprintf(filename,"}\n");

  // Print the matrix
  for (rowLoop in 1:rows){
    clear(rowString);
    rowString=_formatComplex(options,matrix[rowLoop;1]);
    if (columns>1){
      for (columnLoop in 2:columns){
        numberString=_formatComplex(options,matrix[rowLoop;columnLoop]);
        rowString=rowString+columnSeperator+numberString;
        }
      }
    if (rowLoop != rows){
      rowString = rowString+" \\\\\n";
      } else {{
        rowString = rowString+"\n";
        }
      }
    fprintf(filename,rowString);
    }

  // Print the postamble
  fprintf(filename,"\\end{array}\\right]\n");
  fprintf(filename,"\\end{equation}\n");
  return 1;
  };

// The functions _formatComplex() and _formatReal() are used by
// wrileltx() to format the entries of the matrix (the underscore is
// intended to indicate (C-style) that they are internal functions)

// _formatComplex reads a Complex and returns a formatted string.
// We care about sensible formatting so 0+0i, 0-1i, 1+0i 1+0.0032000i
// should be returned as 0, -i, 1 and 1+0.0032i respectively. This makes
// for a messy pile of conditionals but what can you do?  This
// function uses _formatReal() (below) which returns a list containing
// a sign string and an absolute value string.

_formatComplex=function(options,x){

  // local variables
  local(outVal,realVal,imVal,realPart,imPart);

  // format the real and imaginary part of the number
  realPart=_formatReal(options , real(x));
  imPart=_formatReal(options , imag(x));

  if (imPart.val=="0"){
    if (realPart.val=="0"){
      outVal="0";
      } else { {
        if (realPart.sgn=="+"){
          outVal=realPart.val;
          } else { {
            outVal="-"+realPart.val;
            }
          }
        }
      }
    } else { {
      if (realPart.val=="0"){
        if (imPart.val=="1"){
          if (imPart.sgn=="+"){
            outVal=options.imUnit;
            } else { {
              outVal="-"+options.imUnit;
              }
            }
          } else { {
            if (imPart.sgn=="+"){
              outVal=imPart.val+options.imUnit;
              } else { {
                outVal="-"+imPart.val+options.imUnit;
                }
              }
            }
          }
        } else { {
          if (imPart.val=="1"){
            if (realPart.sgn=="+"){
              outVal=realPart.val+imPart.sgn+options.imUnit;
              } else { {
                outVal="-"+realPart.val+imPart.sgn+options.imUnit;
                }
              }
            } else { {
              if (realPart.sgn=="+"){
                outVal=realPart.val+imPart.sgn+imPart.val+options.imUnit;
                } else { {
                  outVal="-"+realPart.val+imPart.sgn+imPart.val+options.imUnit;
                  }
                }
              }
            }
          }
        }
      }
    }
  return outVal;
  };

// _formatReal() returns a list with members sgn ("+" or "-")
// and val (a string of the absolute val of x in a nice format)

_formatReal=function(options,x){

  local(chopped,temp,outSgn,outVal,stend,estart,esign,eastart,exponent,dummy);

  // print out the number to a string temp
  sprintf(temp,options.fstring,x);

  // split it up into its constituent characters
  chopped=strsplt(temp);

  // Where is the end of the string?
  stend=chopped.n;

  // if exponential format then get the exponent
  estart=find(chopped=="e");
  if (estart!=[]){
    // the sign of the exponent
    expsgn=chopped[estart+1];
    // the start of the absolute value of the exponemt
    eastart=estart+2;
    // remove initial zeros
    while(chopped[eastart]=="0" && eastart<chopped.n ){eastart++;}
    // initialize the exponent string
    exponent="";
    for (dummy in eastart:chopped.n){
      exponent=exponent+chopped[dummy];
      }
    stend=estart-1;
    }

  // remove any zeros from the end (but make sure that there is a point first)
  if (any(chopped==".")){
    while (chopped[stend]=="0"){
      stend--;
      }
    }

  // if there is a hanging point remove it
  if (chopped[stend]=="."){ stend--; }

  // write it to the output string

  // find the sign
  if (chopped[1]=="-"){
    outSgn="-";
    ststart=2;
    } else { {
      outSgn="+";
      ststart=1;
      }
    }

  // skip over any blanks
  while (chopped[ststart]==" "){ststart++;}

  // write the absolute value string
  outVal=chopped[ststart];
  if (stend>ststart){
    for (dummy in ststart+1:stend){
      outVal=outVal+chopped[dummy];
      }
    }

  // add the exponent
  if (options.format=="exponential"){
    if (exponent!="0"){
      outVal=outVal+"\\times 10";
      if (expsgn=="+"){
        if (exponent!="1"){
          outVal=outVal+"^{"+exponent+"}";
          }
        } else { {
          outVal=outVal+"^{-"+exponent+"}";
          }
        }
      }
    }
  // ... and end
  return <<val=outVal;sgn=outSgn>>;
  };


