#include <stdio.h> 
#include <time.h> 
#include <pwd.h> 
#include <sys/types.h> 
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "xwin.h"
#include "postscript.h"

#define MAXBUF 256

int layer=0;
int linetype=0;
int pennum=0;
int filltype=0;
int in_line=0;		/* set to 1 when working on a line */
int in_poly=0;

int this_pen=0;		// we need to save state at start of line
			// pen, line and fill

int debug=0;

FILE *fp=NULL;		// file descriptor for output file
OMODE outputtype=POSTSCRIPT;	// 0=postscript, 1=gerber

double bbllx, bblly, bburx, bbury;

void ps_set_outputtype(OMODE mode) {
   if (debug) printf("ps_set_outputtype: %d\n",mode);
   outputtype=mode;
}

void ps_set_file(FILE *fout) 
{
  if (debug) printf("ps_set_file:\n");
  fp=fout; 
}

void ps_set_layer(int layernumber) {
    layer=layernumber;
}

void ps_set_fill(int fill)
{
    filltype=fill;
}

void ps_set_line(int line)
{
    linetype=line;
}

int dxf_penmap[] = {		// input is pen called, output is dxf pen
    0,0,1,3,
    5,4,6,2,
    1,3,5,4,
    6,3,9,8
};

// 0 blk           0
// 1 wht           0
// 2 red           1
// 3 green         3
// 4 blue          5
// 5 aqua          4
// 6 mag           6
// 7 yellow        2
// 8 pink          1
// 9 lime          3
// 10 pale blue    5
// 11 blue-green   4
// 12 purple       6
// 13 khaki        3
// 14 light grey   9
// 15 dark grey    8

void ps_set_pen(int pen) 
{
    if (outputtype==DXF) {
        pennum=dxf_penmap[pen%16];
    } else {
	pennum=pen;
    }
}

/*
Letter  0 0 8.5i 11i 
Legal   0 0 8.5i 14i
Tabloid 0 0  11i 17i
ANSI-C  0 0  17i 22i
ANSI-D  0 0  22i 34i
ANSI-E  0 0  34i 44i
*/


void ps_preamble(dev, prog, pdx, pdy, llx, lly, urx, ury) 
char *dev;
char *prog;
double pdx, pdy; 	    /* page size in inches */
double llx, lly, urx, ury;  /* drawing bounding box in user coords */
{

    time_t timep;
    char buf[MAXBUF];
    double xmid, ymid;
    double xdel, ydel;
    double s1, s2, scale;
    double xmax, ymax;
    int landscape;
    double tmp;
    int debug=0;

    if ( outputtype == AUTOPLOT) {
       fprintf(fp,"back\n");
       fprintf(fp,"isotropic\n");
       return;
    } else if (outputtype == GERBER) {
       fprintf(fp,"G4 Title: %s*\n", dev);
       fprintf(fp,"G4 Creator: %s*\n", prog);
       timep=time(&timep);
       ctime_r(&timep, buf);
       buf[strlen(buf)-1]='\0';
       fprintf(fp,"G4 CreationDate: %s*\n", buf);
       gethostname(buf,MAXBUF);
       fprintf(fp,"G4 For: %s@%s (%s)*\n", getpwuid(getuid())->pw_name, buf, getpwuid(getuid())->pw_gecos );
       fprintf(fp,"G01*\n");	// linear interpolation
       fprintf(fp,"G70*\n");	// inches
       fprintf(fp,"G90*\n");	// absolute
       fprintf(fp,"%%MOIN*%%\n");	// absolute
       fprintf(fp,"G04 Gerber Fmt 3.4, Leading zero omitted, Abs format*\n");
       fprintf(fp,"%%FSLAX34Y34*%%\n");
       fprintf(fp,"G04 Define D10 to be a 5mil circle*\n");
       fprintf(fp,"%%ADD10C,.005*%%\n");
       fprintf(fp,"G04 Make D10 the default aperture*\n");
       fprintf(fp,"G54D10*\n");
       return;
    } else if (outputtype == DXF) {
       fprintf(fp, "0\n");
       fprintf(fp, "SECTION\n");
       fprintf(fp, "2\n");
       fprintf(fp, "ENTITIES\n");
       return;
    }

    bbllx=llx; bblly=lly; bburx=urx; bbury=ury;
    if (debug) printf("llx:%g lly:%g urx:%g ury:%g\n", llx, lly, urx, ury);

    /* now create transform commands based on bb and paper size */
    /* FIXME, check that bb is in canonic order */
    if ((llx > urx) || (lly > ury)) {
    	printf("bounding box error in postscript.c\n");
    }

/*
plotting postcript to zhongsheet.ps
currep min/max bounds: 0, 0, 8.5, 11
currep vp bounds: 0, 0, 8.5, 11
currep screen bounds: 142.932, 787, 751.068, -1.42109e-14
gwidth/height: 894, 787
bounds to ps: 142.932, -1.42109e-14, 751.068, 787
setting landscape
*/
    xmid = (llx+urx)/2.0;
    ymid = (lly+ury)/2.0;
    xdel = 1.05*fabs(urx-llx);
    ydel = 1.05*fabs(ury-lly);

    // pdx=8.5;
    // pdy=11.0;

    xmax=pdx*72.0;
    ymax=pdy*72.0;
    // printf("aspect is xdel: %g, ydel: %g\n", xdel, ydel);
    if (xdel > ydel) { /* flip aspect */
	landscape=1;
    	// printf("setting landscape\n");
	tmp=xmax;  xmax=ymax; ymax=tmp;
    } else {
    	// printf("setting portrait\n");
	landscape=0;
    }

    s1 = xmax/xdel;
    s2 = ymax/ydel;

    scale = s2;
    if (s1 < s2) {
    	scale = s1;
    }

    fprintf(fp,"%%!PS-Adobe-2.0\n");
    fprintf(fp,"%%%%Title: %s\n", dev);
    fprintf(fp,"%%%%Creator: %s\n", prog);
    timep=time(&timep);
    fprintf(fp,"%%%%CreationDate: %s", ctime(&timep));
    gethostname(buf,MAXBUF);
    fprintf(fp,"%%%%For: %s@%s (%s)\n", 
    	getpwuid(getuid())->pw_name, 
	buf,
	getpwuid(getuid())->pw_gecos );
    if (landscape) {
	fprintf(fp,"%%%%Orientation: Landscape\n");
    } else {
	fprintf(fp,"%%%%Orientation: Portrait\n");
    }
    fprintf(fp,"%%%%Pages: (atend)\n");
    fprintf(fp,"%%%%BoundingBox: 0 0 %d %d\n", 
    	(int) (pdx*72.0), (int) (pdy*72.0));
    fprintf(fp,"%%%%DocumentPaperSizes: custom\n");
    fprintf(fp,"%%%%Magnification: 1.0000\n");
    fprintf(fp,"%%%%EndComments\n");
    fprintf(fp,"%%%%BeginSetup\n");
    fprintf(fp,"mark {\n");
    fprintf(fp,"%%BeginFeature: *PageRegion custom\n");
    fprintf(fp,"<</PageSize [ %d %d ] /ImagingBBox null >> setpagedevice\n",
    	(int)( pdx*72.0), (int) (pdy*72.0));
    fprintf(fp,"%%EndFeature\n");
    fprintf(fp,"} stopped cleartomark\n");
    fprintf(fp,"%%%%EndSetup\n");
    fprintf(fp,"%%%%BeginProlog\n");

    fprintf(fp,"/$Pig2psDict 200 dict def\n");
    fprintf(fp,"$Pig2psDict begin\n");
    fprintf(fp,"$Pig2psDict /mtrx matrix put\n");
    fprintf(fp,"/c-1 {0 setgray} bind def\n");
    fprintf(fp,"/c0 {0.000 0.000 0.000 srgb} bind def\n");  	// blk
    fprintf(fp,"/c1 {0.000 0.000 0.000 srgb} bind def\n");  	// blk
    fprintf(fp,"/c2 {1.000 0.000 0.000 srgb} bind def\n");	// red
    fprintf(fp,"/c3 {0.000 1.000 0.000 srgb} bind def\n");	// grn
    fprintf(fp,"/c4 {0.000 0.000 1.000 srgb} bind def\n");	// blu
    fprintf(fp,"/c5 {0.000 1.000 1.000 srgb} bind def\n");	// aqu
    fprintf(fp,"/c6 {1.000 0.000 1.000 srgb} bind def\n");	// mag
    fprintf(fp,"/c7 {1.000 1.000 0.000 srgb} bind def\n");	// yel
    fprintf(fp,"/c8 {1.000 0.500 0.500 srgb} bind def\n");	// pnk
    fprintf(fp,"/c9 {0.500 1.000 0.500 srgb} bind def\n");	// lime
    fprintf(fp,"/c10 {0.500 0.500 1.000 srgb} bind def\n");	// pale blu
    fprintf(fp,"/c11 {0.000 0.500 0.500 srgb} bind def\n");	// blu-grn
    fprintf(fp,"/c12 {0.500 0.000 0.500 srgb} bind def\n");	// pur
    fprintf(fp,"/c13 {0.500 0.500 0.000 srgb} bind def\n");	// khaki
    fprintf(fp,"/c14 {0.700 0.700 0.700 srgb} bind def\n");	// lt gr
    fprintf(fp,"/c15 {0.400 0.400 0.400 srgb} bind def\n");	// dk gr
    fprintf(fp,"/c16 {1.000 0.000 1.000 srgb} bind def\n");
    fprintf(fp,"/c17 {1.000 1.000 0.000 srgb} bind def\n");
    fprintf(fp,"/c18 {0.000 0.000 0.000 srgb} bind def\n");
    fprintf(fp,"/c19 {0.300 0.300 0.300 srgb} bind def\n");
    fprintf(fp,"/c20 {0.600 0.600 0.600 srgb} bind def\n");
    fprintf(fp,"/cp {closepath} bind def\n");
    fprintf(fp,"/ef {eofill} bind def\n");
    fprintf(fp,"/gr {grestore} bind def\n");
    fprintf(fp,"/gs {gsave} bind def\n");
    fprintf(fp,"/sa {save} bind def\n");
    fprintf(fp,"/rs {restore} bind def\n");
    fprintf(fp,"/kc { currentrgbcolor [/Pattern /DeviceRGB] setcolorspace } bind def\n");
    fprintf(fp,"/l {lineto} bind def\n");
    fprintf(fp,"/m {moveto} bind def\n");
    fprintf(fp,"/rm {rmoveto} bind def\n");
    fprintf(fp,"/n {newpath} bind def\n");
    fprintf(fp,"/s {stroke} bind def\n");
    fprintf(fp,"/sh {show} bind def\n");
    fprintf(fp,"/slc {setlinecap} bind def\n");
    fprintf(fp,"/slj {setlinejoin} bind def\n");
    fprintf(fp,"/slw {setlinewidth} bind def\n");
    fprintf(fp,"/srgb {setrgbcolor} bind def\n");
    fprintf(fp,"/rot {rotate} bind def\n");
    fprintf(fp,"/sc {scale} bind def\n");
    fprintf(fp,"/sd {setdash} bind def\n");
    fprintf(fp,"/ff {findfont} bind def\n");
    fprintf(fp,"/sf {setfont} bind def\n");
    fprintf(fp,"/scf {scalefont} bind def\n");
    fprintf(fp,"/sw {stringwidth} bind def\n");
    fprintf(fp,"/tr {translate} bind def\n");
    fprintf(fp,"end\n");
    fprintf(fp,"save\n");

    // fprintf(fp,"newpath\n");
    // fprintf(fp,"0 %d moveto\n", (int) (pdy*72.0));
    // fprintf(fp,"0 0 lineto\n");
    // fprintf(fp,"%d 0 lineto\n", (int) (pdx*72.0));
    // fprintf(fp,"%d %d lineto closepath clip newpath\n", 
    // 	(int) (pdx*72.0), (int) (pdy*72.0));
    // /* fprintf(fp,"%%39.9 49.0 translate\n"); */

    fprintf(fp,"/$Pig2psBegin\n");
    fprintf(fp,"{$Pig2psDict begin /$Pig2psEnteredState save def}def\n");
    fprintf(fp,"/$Pig2psEnd {$Pig2psEnteredState restore end} def\n");

    fprintf(fp,"%%BeginPageSetup\n");
    fprintf(fp,"%%BB is %g,%g %g,%g\n", llx, lly, urx, ury);	
    if (landscape) {
    	fprintf(fp,"%g %g scale\n", scale, scale);
	fprintf(fp,"%g %g translate\n", 
	    (ymax/(2.0*scale))+ymid, (xmax/(2.0*scale))-xmid);
	fprintf(fp,"90 rotate\n");
    } else {
    	fprintf(fp,"%g %g scale\n", scale, scale);
	fprintf(fp,"%g %g translate\n",
	    (xmax/(2.0*scale))-xmid,  (ymax/(2.0*scale))-ymid);
    }
    fprintf(fp,"%%EndPageSetup\n");
    fprintf(fp,"%%%%EndProlog\n");
    fprintf(fp,"$Pig2psBegin\n");
    fprintf(fp,"10 setmiterlimit\n");
    fprintf(fp,"1 slj 1 slc\n");
    fprintf(fp,"1.0 slw\n");

//
//    fprintf(fp,"%% start stipple patterns\n");
//    int i,j; char *ps;
//    for (i=0;(ps=get_stipple_bits(i))!=NULL;i++) {
//        fprintf(fp,"/p%d {\n",i);
//        fprintf(fp,"8 8 true  matrix {<");
//	for (j=0; j<=7; j++) {
//    	    fprintf(fp,"%02x",ps[j]&255);
//	}
//        fprintf(fp,">} imagemask\n");
//    	fprintf(fp,"} bind def\n");
//	fprintf(fp,"<< /PatternType 1\n");
//        fprintf(fp,"   /PaintType 2\n");
//	fprintf(fp,"   /TilingType 1\n");
//	fprintf(fp,"   /XStep 8 /YStep 8\n");
//	fprintf(fp,"   /BBox [0 0 8 8]\n");
//	fprintf(fp,"   /PaintProc { p%d }\n",i);
//	fprintf(fp,">>\n");
//	fprintf(fp,"matrix makepattern /f%d exch def\n",i);
//   }
//  fprintf(fp,"%% end stipple patterns\n");


    fprintf(fp,"%% here starts figure;\n");
}


void ps_end_line()
{
    extern int this_pen;
    extern int in_line;

    if (debug) printf("ps_end_line:\n");

    if (outputtype == POSTSCRIPT) {
	fprintf(fp, "c%d\n", this_pen);
	fprintf(fp, "s\n");
    }
    in_poly=0;
    in_line=0;
}

double xold, yold;
int in_progress=0;

void ps_start_line(x1, y1)
double x1, y1;
{
    extern int pennum;
    extern int linetype;
    extern int in_line;
    extern int this_pen;

    if (debug) printf("ps_start_line:\n");

    if (in_line) {
    	ps_end_line(fp);
    } 
    in_line++;

    this_pen=pennum;		/* save characteristics at start of line */

    if (outputtype == POSTSCRIPT) {
	fprintf(fp, "n\n");
        // flip y coordinate from Xwin coords
        // y1 = bblly-(y1-bbury);
	fprintf(fp, "%g %g m\n",x1, y1);
    } else if (outputtype == AUTOPLOT) {		// AUTOPLOT
	fprintf(fp, "jump\n");
	fprintf(fp, "%g %g\n",x1, y1);
    } else if (outputtype == GERBER) {			// GERBER
        if (!in_poly) {
	   fprintf(fp,"G04 LAYER %d *\n", layer);
	}
	fprintf(fp, "X%04dY%04dD02*\n",(int)((x1*10000.0)+0.5), (int)((y1*10000.0)+0.5));
    } else if (outputtype == DXF) {
	fprintf(fp, "0\n");		// start of entity
	fprintf(fp, "LINE\n");		// line type entity
	fprintf(fp, "8\n");		// layer name to follow
	fprintf(fp, "%d\n",layer);	// layer name 
	fprintf(fp, "62\n");		// color flag
	fprintf(fp, "%d\n", pennum);	// color
	fprintf(fp, "10\n%g\n", x1);	// initial x value
	fprintf(fp, "20\n%g\n", y1);	// initial y value
	in_progress=0;
    }
}


void ps_continue_line(x1, y1)
double x1, y1;
{
    if (debug) printf("ps_continue_line:\n");

    if (outputtype == POSTSCRIPT) {
	// flip y coordinate from Xwin coords
	// y1 = bblly-(y1-bbury);
	fprintf(fp, "%g %g l\n",x1, y1);
    } else if (outputtype == AUTOPLOT) {		// AUTOPLOT
	fprintf(fp, "%g %g\n",x1, y1);
    } else if (outputtype == GERBER) {			// GERBER
	fprintf(fp, "X%04dY%04dD01*\n",(int)((x1*10000.0)+0.5), (int)((y1*10000.0)+0.5));
    } else if (outputtype == DXF) {
	if (!in_progress) {
	    fprintf(fp, "11\n%g\n", x1);	// initial x value
	    fprintf(fp, "21\n%g\n", y1);	// initial y value
	    in_progress++;
	} else {
	    fprintf(fp, "0\n");		// start of entity
	    fprintf(fp, "LINE\n");	// line type entity
	    fprintf(fp, "8\n");		// layer name to follow
	    fprintf(fp, "%d\n",layer);	// layer name 
	    fprintf(fp, "62\n");		// color flag
	    fprintf(fp, "%d\n", pennum);	// color
	    fprintf(fp, "10\n%g\n", xold);	// initial x value
	    fprintf(fp, "20\n%g\n", yold);	// initial y value
	    fprintf(fp, "11\n%g\n", x1);	// initial x value
	    fprintf(fp, "21\n%g\n", y1);	// initial y value
	}
	xold = x1;
	yold = y1;
    }
    in_line++;
}

void ps_start_poly() 
{
    if (debug) printf("ps_start_poly:\n");
    if (outputtype == GERBER) {			// GERBER
       fprintf(fp,"G04 LAYER %d *\n", layer);
       fprintf(fp, "G36*\n");
    }
    in_poly++;
}

void ps_end_poly() 
{
    if (debug) printf("ps_end_poly:%d\n",outputtype);
    if (outputtype == GERBER) {			// GERBER
       fprintf(fp, "G37*\n");
    }
    in_poly=0;
}

void ps_postamble()
{
    if (in_line) {
    	ps_end_line(fp);
	in_line=0;
    } 
    if (outputtype == POSTSCRIPT) {
	fprintf(fp,"%% here ends figure;\n");
	fprintf(fp,"%%Pig2psEnd\n");
	fprintf(fp,"rs\n");
	fprintf(fp,"showpage\n");
	fprintf(fp,"%%%%EOF\n");
    } else if (outputtype == GERBER) {
        fprintf(fp,"G04 LAYER 999 *\n");
	fprintf(fp,"M02*\n");
    } else if (outputtype == DXF) {
        fprintf(fp,"0\nENDSEC\n");
        fprintf(fp,"0\nEOF\n");
    }
    fclose(fp);
}
