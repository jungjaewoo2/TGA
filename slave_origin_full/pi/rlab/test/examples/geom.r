//-------------------------------------------------------------------//
// Synopsis: Rlab -> Geomview Interface

// Syntax:   geommesh ( Mesh-List , Color )
//           geomclose ( )

// Description:

//     These functions provide a crude, but effective Geomview
//     interface to Rlab. This is work in progress.

//     geommesh takes two arguments:
//        Mesh-List is a list containing (at least) 3 elements:
//        x       The x-coordinates
//        y       The y-coordinates
//        z       The mesh elevations.
//     Color, an optional argument, is non-zero (default) to instruct
//     geommesh to code the elevations with a color-scale. If color is
//     zero, then the mesh is drawn in greyscale.

// Ian Searle, 11/11/95
//-------------------------------------------------------------------//

#
# The Geomview process we will pipe data to.
#

static (GEOM)
GEOM = "|/usr/local/bin/geomview -c -";

#
# The object tag.
#

static (objt)
objt = 1;

# geommesh ( Mesh-List )
#
# ML:  The input Mesh-List. ML must have (at a minimum) elements:
#      x   The x-coordinates
#      y   The y-coordinates
#      z   The elevations
#

geommesh = function ( ML, color )
{
  if (!exist (color)) { color = 1; }

  #
  # Write out the header info...
  #

  fprintf (GEOM, "(geometry rlab-mesh-%i { : rlab-mesh-%i })\n", objt, objt);
  fprintf (GEOM, "(read geometry { define rlab-mesh-%i\n", objt);
  objt++;

  if (color)
  {
    fprintf (GEOM, "CMESH\n");
  else
    fprintf (GEOM, "MESH\n");
  }
  fprintf (GEOM, "%i  %i\n", ML.x.n, ML.y.n);

  #
  # Figure out the color scales...
  #

  MAX = max (max (ML.z));
  MIN = min (min (ML.z));
  range = MAX - MIN;

  #
  # Now start writing out the vertices...
  #

  for (j in 1:ML.y.n)
  {
    for (i in 1:ML.x.n)
    {
      # Write out the vertex values...
      fprintf (GEOM, "%f %f %f  ", ML.x[i], ML.y[j], ML.z[i;j]);

      if (color)
      {
	# Write out the color values...
	val = (ML.z[i;j] - MIN) / range;
	rcolor = 1.58*(1 - exp (-val));
	gcolor = val.^2;
	bcolor = 0;
	fprintf (GEOM, "%f %f %f 1.0\n", rcolor, gcolor, bcolor);
      else
	fprintf (GEOM, "\n");
      }
    }
    fprintf (GEOM, "\n");
  }

  fprintf (GEOM, "})\n");
};

#
# Close the Geomview process.
#

geomclose = function ( )
{
  # Reset the object tag
  objt = 1;
  return close (GEOM);
};
