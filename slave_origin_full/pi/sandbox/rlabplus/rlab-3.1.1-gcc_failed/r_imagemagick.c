//
// rimagemagick.c:  rlab's interface to ImageMagick
// Marijan Kostrun, VII-VIII-2013
//
// This file is a part of RLaB + rlabplus
// Copyright (C) 2013  Marijan Kostrun
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// See the file ../COPYING

#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "symbol.h"
#include "mem.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mds.h"
#include "list.h"
#include "btree.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"

#include <stdio.h>
#include <fcntl.h>

//
#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"


#ifdef HAVE_IMAGEMAGICK  /* Set in config.h */

#define THIS_FILE "r_imagemagick.c"

#include <wand/MagickWand.h>

#define MAX_NUMBER_MAGICK_WANDS 32
static int InitMagickWandLib=0;
static int default_wand_idx=-1;
static MagickWand *magick_wand[MAX_NUMBER_MAGICK_WANDS] = {0};
static unsigned char im_debug=1;
static ExceptionType severity;
static char *description=0;

#define RLAB_MAGICK_QUANTUM_SIZE 255

#define MAX_LEN_DEFAULT_TEXT 1024
static char default_text[MAX_LEN_DEFAULT_TEXT];

static void print_recognized_methods(void)
{
  fprintf(stderr, "Implemented methods are: "
  "adaptiveblur,"
  "adaptiveresize,"
  "sharpen,"
  "adaptivesharpen,"
  "adaptivethreshold,"
  "blackthreshold,"
  "threshold,"
  "addnoise,"
  "affine (not working),"
  "annotate,"
  "autogamma,"
  "autolevel,"
  "blueshift,"
  "gaussianblur,"
  "motionblur,"
  "radialblur,"
  "selectiveblur,"
  "blur,"
  "border,"
  "brightnesscontrast,"
  "charcoal,"
  "chop,"
  "clamp,"
  "clip,"
  "cdl,"
  "colormatrix,"
  "colorimage,"
  "colorize,"
  "contraststretch,"
  "sigmoidalcontrast,"
  "contrast,"
  "convolve,"
  "crop,"
  "cyclecolormap,"
  "deskew,"
  "despecle,"
  "edge,"
  "emboss,"
  "enhance,"
  "equalize,"
  "extent,"
  "flip,"
  "flop,"
  "gamma,"
  "implode,"
  "level,"
  "linearstretch,"
  "liquidrescale,"
  "modulate,"
  "negate,"
  "normalize,"
  "oilpaint,"
  "polaroid,"
  "posterize,"
  "randomthreshold,"
  "resize,"
  "roll,"
  "rotate,"
  "sample,"
  "scale,"
  "sepia,"
  "shade,"
  "shadow,"
  "shave,"
  "shear,"
  "sketch,"
  "solarize,"
  "splice,"
  "spread,"
  "strip,"
  "swirl,"
  "thumbnail,"
  "transparentpaint,"
  "transpose,"
  "transverse,"
  "trim,"
  "uniquecolors,"
  "unsharpmask,"
  "vignette,"
  "wave,"
  "whitethreshold"
  "\n");
}

static void WriteWandException(MagickWand *wand)
{
  char *d;

  d = MagickGetException(wand,&severity);

  if (description)
    GC_free(description);

  description = cpstr(d);
  if (im_debug)
    fprintf(stderr,"%s %s %lu %s\n",GetMagickModule(),d);

  d=(char *) MagickRelinquishMemory(d);

  return;
}

//
// wrapture is here ! repent!
//    Or maybe not. Just wrappers for the functions so that they accept double arguments:
//
static MagickBooleanType _MagickAdaptiveResizeImage(MagickWand *wand,double columns, double rows)
{
  size_t r,c;
  if (columns <=1)
    c = (size_t) ((double) MagickGetImageWidth(wand) * columns);
  else
    c = (size_t) columns;
  if (rows <= 1)
    r = (size_t) ((double) MagickGetImageHeight(wand) * rows);
  else
    r = (size_t) rows;
  return MagickAdaptiveResizeImage((MagickWand *)wand, c, r);
}

static MagickBooleanType _MagickAdaptiveThresholdImage(MagickWand *wand,
                                                       double width,double height,double offset)
{
  return MagickAdaptiveThresholdImage((MagickWand *)wand,
                              (const size_t) width,
                              (const size_t) height,
                              (const size_t) offset);
}
static MagickBooleanType _MagickCycleColormapImage(MagickWand *wand, double displace)
{
  return MagickCycleColormapImage((MagickWand *)wand, (const ssize_t) displace);
}
static MagickBooleanType _MagickContrastImage(MagickWand *wand, double sharpen)
{
  return MagickContrastImage((MagickWand *)wand, (const MagickBooleanType) sharpen);
}
static MagickBooleanType _MagickNegateImage(MagickWand *wand, double gray)
{
 return MagickNegateImage((MagickWand *)wand, (const MagickBooleanType) gray);
};
static MagickBooleanType _MagickNegateImageChannel(MagickWand *wand, const ChannelType channel,
                                                   double gray)
{
  return MagickNegateImageChannel((MagickWand *)wand, (const ChannelType) channel,
                                  (const MagickBooleanType) gray);
}
static MagickBooleanType _MagickPolaroidImage(MagickWand *wand,double angle)
{
  MagickBooleanType rval;
  DrawingWand *drawing_wand = NewDrawingWand();
  rval = MagickPolaroidImage((MagickWand *)wand,
                             (const DrawingWand *) drawing_wand,
                             (const double) angle);
  drawing_wand = DestroyDrawingWand(drawing_wand);
  return rval;
}
static MagickBooleanType _MagickPosterizeImage(MagickWand *wand,double levels,double dither)
{
  return MagickPosterizeImage((MagickWand *)wand, (const size_t) levels, (const MagickBooleanType) dither);
}
static MagickBooleanType _MagickRollImage(MagickWand *wand, double x, double y)
{
  return MagickRollImage((MagickWand *)wand, (const size_t) x, (const size_t) y);
}
static MagickBooleanType _MagickSampleImage(MagickWand *wand, double columns, double rows)
{
  return MagickSampleImage((MagickWand *)wand, (const size_t) columns, (const size_t) rows);
}
static MagickBooleanType _MagickScaleImage(MagickWand *wand, double columns, double rows)
{
  return MagickScaleImage((MagickWand *)wand, (const size_t) columns, (const size_t) rows);
}
static MagickBooleanType _MagickShadeImage(MagickWand *wand, double gray, double azimuth, double elevation)
{
  return MagickShadeImage((MagickWand *)wand, (const MagickBooleanType) gray,
                          (const double) azimuth, (const double) elevation);
}
static MagickBooleanType _MagickShadowImage(MagickWand *wand, double opacity,
                                            double sigma, double x, double y)
{
  return MagickShadowImage((MagickWand *)wand, (const double) opacity,
                           (const double)sigma, (const size_t) x,
                           (const size_t) y);
}
static MagickBooleanType _MagickShaveImage(MagickWand *wand, double columns, double rows)
{
  return MagickShaveImage((MagickWand *)wand, (const size_t) columns, (const size_t) rows);
}
static MagickBooleanType _MagickSigmoidalContrastImage(MagickWand *wand, double sharpen,
                                                       double alpha, double beta)
{
  return MagickSigmoidalContrastImage((MagickWand *) wand,
                                      (const MagickBooleanType) sharpen,
                                      (const double) alpha,
                                      (const double) beta);
}
static MagickBooleanType _MagickSigmoidalContrastImageChannel(MagickWand *wand, ChannelType channel,
                                                              double sharpen, double alpha, double beta)
{
  return MagickSigmoidalContrastImageChannel((MagickWand *) wand, (const ChannelType) channel,
                                      (const MagickBooleanType) sharpen,
                                      (const double) alpha,
                                      (const double) beta);
}
static MagickBooleanType _MagickSpliceImage(MagickWand *wand,
                                            double width, double height, double x, double y)
{
  return MagickSpliceImage((MagickWand *)wand,
                           (const size_t) width, (const size_t) height, (const size_t) x,
                           (const size_t) y);
}
static MagickBooleanType _MagickThumbnailImage(MagickWand *wand, double columns, double rows)
{
  return MagickThumbnailImage((MagickWand *)wand,(const size_t) columns,(const size_t) rows);
}
static MagickBooleanType _MagickVignetteImage(MagickWand *wand, double black_point, double white_point,
                                              double x,double y)
{
  return MagickVignetteImage((MagickWand *)wand,
                      (const double) black_point,
                      (const double) white_point,
                      (const size_t) x,
                      (const size_t) y);
}



//
// set default wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".method"
#define IM_MAX_NO_METHOD_PARAMS 12
Ent *
ent_im_wand_apply (int nargs, Datum args[])
{
  Ent *e1=0, *en=0, *rent;
  int i,j,rval=0;
  char *method=0;
  int np=0, have_ch=0;  // for counting parameters to the methods
  char *lb,*rb, *c=0;
  unsigned char ar=0,ac=0,ag=0,am=0,ab=0,ay=0,ao=0,aB=0,aM=0;
  double params[IM_MAX_NO_METHOD_PARAMS] = {0};
  MDR *mdr_param=0;

  // parameters needed for calling the functions
  MagickBooleanType (*func)() = 0;      // function name
  MagickBooleanType (*func_ch)() = 0;   // function name
  NoiseType noise_type;                 // special function parameters

  MagickBooleanType bval;

  if (nargs < 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Apply method to active image in default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(method,param1,....),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'method' is the name of the method to be applied, while.\n");
    fprintf (stdout,
             THIS_SOLVER ":   'param1,..' are the method parameters.\n");
    print_recognized_methods();
    fprintf (stdout,
             THIS_SOLVER ": Check manual for method parameters.\n");
    rerror ( THIS_SOLVER ": requires 1 or more arguments");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs >= 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_STRING)
      method  = class_char_pointer(e1);
  }

  if (!method)
    rerror(THIS_SOLVER  ": First argument 'method' is required!\n");

  // identify the method and its arguments
  if (strstr(method,"adaptiveblur"))
  {
    // parameters: radius, sigma
    have_ch = 1;
    func = &MagickAdaptiveBlurImage;
    func_ch = &MagickAdaptiveBlurImageChannel;
    np = 2;
  }
  else if (strstr(method,"adaptiveresize"))
  {
    // parameters: width, height or fractional (<=1) of original size
    have_ch = 0;
    func = &_MagickAdaptiveResizeImage;
    np = 2;
  }
  else if (strstr(method,"sharpen"))
  {
    // parameters: radius, sigma
    have_ch = 1;
    func = &MagickSharpenImage;
    func_ch = &MagickSharpenImageChannel;
    np = 2;
  }
  else if (strstr(method,"adaptivesharpen"))
  {
    // parameters: radius, sigma
    have_ch = 1;
    func = &MagickAdaptiveSharpenImage;
    func_ch = &MagickAdaptiveSharpenImageChannel;
    np = 2;
  }
  else if (strstr(method,"adaptivethreshold"))
  {
    // parameters: width, height, offset
    have_ch = 0;
    func = &_MagickAdaptiveThresholdImage;
    np = 3;
  }
  else if (strstr(method,"blackthreshold"))
  {
    // parameters: r,g,b,o
    have_ch = 0;
    func = &MagickBlackThresholdImage;
    np = -1;
  }
  else if (strstr(method,"threshold"))
  {
    // parameters: threshold
    have_ch = 1;
    func = &MagickThresholdImage;
    func_ch = &MagickThresholdImageChannel;
    np = 1;
  }
  else if (strstr(method,"addnoise"))
  {
    // parameters: "type", attanuate
    have_ch = 1;
    func = &MagickAddNoiseImage;
    func_ch = &MagickAddNoiseImageChannel;
    np = 2;
  }
  else if (strstr(method,"affine"))
  {
    // parameters: rx,sx,sy,ry,tx,ty
    have_ch = 0;
    func = &MagickAffineTransformImage;
    np = 6;
  }
  else if (strstr(method,"annotate"))
  {
    // parameters: text, [x,y], angle
    have_ch = 0;
    func = &MagickAnnotateImage;
    np = -1;
  }
  else if (strstr(method,"autogamma"))
  {
    have_ch = 1;
    func = &MagickAutoGammaImage;
    func_ch = &MagickAutoGammaImageChannel;
    np = 0;
  }
  else if (strstr(method,"autolevel"))
  {
    have_ch = 1;
    func = &MagickAutoLevelImage;
    func_ch = &MagickAutoLevelImageChannel;
    np = 0;
  }
  else if (strstr(method,"blueshift"))
  {
    // parameters: factor
    have_ch = 0;
    func = &MagickBlueShiftImage;
    np = 1;
  }
  else if (strstr(method,"gaussianblur"))
  {
    // parameters: radius,sigma
    have_ch = 1;
    func = &MagickGaussianBlurImage;
    func_ch = &MagickGaussianBlurImageChannel;
    np = 2;
  }
  else if (strstr(method,"motionblur"))
  {
    // parameters: radius, sigma, angle
    have_ch = 1;
    func = &MagickMotionBlurImage;
    func_ch = &MagickMotionBlurImageChannel;
    np = 3;
  }
  else if (strstr(method,"radialblur"))
  {
    // parameters: angle
    have_ch = 1;
    func = &MagickRadialBlurImage;
    func_ch = &MagickRadialBlurImageChannel;
    np = 1;
  }
  else if (strstr(method,"selectiveblur"))
  {
    // parameters: radius,sigma,threshold
    have_ch = 1;
    func = &MagickSelectiveBlurImage;
    func_ch = &MagickSelectiveBlurImageChannel;
    np = 3;
  }
  else if (strstr(method,"blur"))
  {
    // parameters: radius, sigma
    have_ch = 1;
    func = &MagickBlurImage;
    func_ch = &MagickBlurImageChannel;
    np = 2;
  }
  else if (strstr(method,"border"))
  {
    // parameters: r,g,b,o,w,h
    have_ch = 0;
    func = &MagickBorderImage;
    np = 6;
  }
  else if (strstr(method,"brightnesscontrast"))
  {
    // parameters: brightness, contrast
    have_ch = 1;
    func = &MagickBrightnessContrastImage;
    func_ch = &MagickBrightnessContrastImageChannel;
    np = 2;
  }
  else if (strstr(method,"charcoal"))
  {
    // parameters: radius, sigma
    have_ch = 0;
    func = &MagickCharcoalImage;
    np = 2;
  }
  else if (strstr(method,"chop"))
  {
    // parameters: width, height, x-offset, y-offset
    have_ch = 0;
    func = &MagickChopImage;
    np = 4;
  }
  else if (strstr(method,"cdl"))
  {
    // parameters: filename, or
    // offset[1x3], slope[1x3], power[1x3], saturation
    have_ch = 0;
    func = &MagickColorDecisionListImage;
    np = -1;
  }
  else if (strstr(method,"colorize"))
  {
    // parameters: filename, or
    // offset[1x3], slope[1x3], power[1x3], saturation
    have_ch = 0;
    func = &MagickColorizeImage;
    np = -1;
  }
  else if (strstr(method,"colormatrix"))
  {
    // parameters: filename, or
    // offset[1x3], slope[1x3], power[1x3], saturation
    have_ch = 0;
    func = &MagickColorMatrixImage;
    np = -1;
  }
  else if (strstr(method,"clamp"))
  {
    have_ch = 1;
    func = &MagickClampImage;
    func_ch = &MagickClampImageChannel;
    np = 0;
  }
  else if (strstr(method,"clip"))
  {
    have_ch = 0;
    func    = &MagickClipImage;
    np = 0;
  }
  else if (strstr(method,"contraststretch"))
  {
    // parameters: black,white
    have_ch = 1;
    func    = &MagickContrastStretchImage;
    func_ch = &MagickContrastStretchImageChannel;
    np = 2;
  }
  else if (strstr(method,"sigmoidalcontrast"))
  {
    // parameters: sharpen, alpha, beta
    have_ch = 1;
    func    = &_MagickSigmoidalContrastImage;
    func_ch = &_MagickSigmoidalContrastImageChannel;
    np = 3;
  }
  else if (strstr(method,"convolve"))
  {
    have_ch = 1;
    func    = &MagickConvolveImage;
    func    = &MagickConvolveImageChannel;
    np = -1;
  }
  else if (strstr(method,"contrast"))
  {
    have_ch = 0;
    func    = &_MagickContrastImage;
    np = 1;
  }
  else if (strstr(method,"crop"))
  {
    // parameters: width, height, x-offset, y-offset
    have_ch = 0;
    func = &MagickCropImage;
    np = 4;
  }
  else if (strstr(method,"cyclecolormap"))
  {
    // parameters: d
    have_ch = 0;
    func = &_MagickCycleColormapImage;
    np = 1;
  }
  else if (strstr(method,"deskew"))
  {
    // parameters: threshold
    have_ch = 0;
    func = &MagickDeskewImage;
    np = 1;
  }
  else if (strstr(method,"despecle"))
  {
    // parameters: threshold
    have_ch = 0;
    func = &MagickDespeckleImage;
    np = 0;
  }
  else if (strstr(method,"edge"))
  {
    // parameters: threshold
    have_ch = 0;
    func = &MagickEdgeImage;
    np = 1;
  }
  else if (strstr(method,"emboss"))
  {
    // parameters: radius, sigma
    have_ch = 0;
    func = &MagickEmbossImage;
    np = 2;
  }
  else if (strstr(method,"enhance"))
  {
    have_ch = 0;
    func = &MagickEnhanceImage;
    np = 0;
  }
  else if (strstr(method,"equalize"))
  {
    have_ch = 1;
    func = &MagickEqualizeImage;
    func_ch = &MagickEqualizeImageChannel;
    np = 0;
  }
  else if (strstr(method,"extent"))
  {
    // parameters: w,h,x,y
    have_ch = 0;
    func = &MagickExtentImage;
    np = 4;
  }
  else if (strstr(method,"flip"))
  {
    have_ch = 0;
    func = &MagickFlipImage;
    np = 0;
  }
  else if (strstr(method,"flop"))
  {
    have_ch = 0;
    func = &MagickFlopImage;
    np = 0;
  }
  else if (strstr(method,"gamma"))
  {
    // parameters: gamma
    have_ch = 1;
    func = &MagickGammaImage;
    func_ch = &MagickGammaImageChannel;
    np = 1;
  }
  else if (strstr(method,"implode"))
  {
    // parameters: gamma
    have_ch = 0;
    func = &MagickImplodeImage;
    np = 1;
  }
  else if (strstr(method,"level"))
  {
    // parameters: gamma
    have_ch = 1;
    func = &MagickLevelImage;
    func_ch = &MagickLevelImageChannel;
    np = 3;
  }
  else if (strstr(method,"linearstretch"))
  {
    // parameters: gamma
    have_ch = 0;
    func = &MagickLinearStretchImage;
    np = 2;
  }
  else if (strstr(method,"liquidrescale"))
  {
    // parameters: cols, rows, delta_x, rigidity
    have_ch = 0;
    func = &MagickLiquidRescaleImage;
    np = 4;
  }
  else if (strstr(method,"modulate"))
  {
    // parameters: brightness, saturation, hue
    have_ch = 0;
    func = &MagickModulateImage;
    np = 3;
  }
  else if (strstr(method,"negate"))
  {
    // parameters: gray
    have_ch = 1;
    func = &_MagickNegateImage;
    func_ch = &_MagickNegateImageChannel;
    np = 1;
  }
  else if (strstr(method,"normalize"))
  {
    have_ch = 1;
    func = &MagickNormalizeImage;
    func_ch = &MagickNormalizeImageChannel;
    np = 0;
  }
  else if (strstr(method,"oilpaint"))
  {
    // parameters: radius
    have_ch = 0;
    func = &MagickOilPaintImage;
    np = 1;
  }
  else if (strstr(method,"polaroid"))
  {
    // parameters: radius
    have_ch = 0;
    func = &_MagickPolaroidImage;
    np = 1;
  }
  else if (strstr(method,"posterize"))
  {
    // parameters: levels, dither
    have_ch = 0;
    func = &_MagickPosterizeImage;
    np = 2;
  }
  else if (strstr(method,"randomthreshold"))
  {
    // parameters: low, high
    have_ch = 1;
    func = &MagickRandomThresholdImage;
    func_ch = &MagickRandomThresholdImageChannel;
    np = 2;
  }
  else if (strstr(method,"resize"))
  {
    // parameters: low, high
    have_ch = 0;
    func = &MagickResizeImage;
    np = 4;
  }
  else if (strstr(method,"roll"))
  {
    // parameters: low, high
    have_ch = 0;
    func = &_MagickRollImage;
    np = 2;
  }
  else if (strstr(method,"rotate"))
  {
    // parameters: low, high
    have_ch = 0;
    func = &MagickRotateImage;
    np = -1;
  }
  else if (strstr(method,"sample"))
  {
    // parameters: low, high
    have_ch = 0;
    func = &_MagickSampleImage;
    np = 2;
  }
  else if (strstr(method,"scale"))
  {
    // parameters: low, high
    have_ch = 0;
    func = &_MagickScaleImage;
    np = 2;
  }
  else if (strstr(method,"sepia"))
  {
    // parameters: threshold
    have_ch = 0;
    func = &MagickSepiaToneImage;
    np = 1;
  }
  else if (strstr(method,"shade"))
  {
    // parameters: gray, azimuth, elevation
    have_ch = 0;
    func = &_MagickShadeImage;
    np = 3;
  }
  else if (strstr(method,"shadow"))
  {
    // parameters: opacity, sigma, x, y
    have_ch = 0;
    func = &_MagickShadowImage;
    np = 4;
  }
  else if (strstr(method,"shave"))
  {
    // parameters: opacity, sigma, x, y
    have_ch = 0;
    func = &_MagickShaveImage;
    np = 2;
  }
  else if (strstr(method,"shear"))
  {
    // parameters: opacity, sigma, x, y
    have_ch = 0;
    func = &MagickShearImage;
    np = -1;
  }
  else if (strstr(method,"sketch"))
  {
    // parameters: radius, sigma, angle
    have_ch = 0;
    func = &MagickSketchImage;
    np = 3;
  }
  else if (strstr(method,"solarize"))
  {
    // parameters: threshold
    have_ch = 0;
    func    = &MagickSolarizeImage;
    np = 1;
  }
  else if (strstr(method,"splice"))
  {
    // parameters: width, height, x, y
    have_ch = 0;
    func    = &_MagickSpliceImage;
    np = 4;
  }
  else if (strstr(method,"spread"))
  {
    // parameters: radius
    have_ch = 0;
    func    = &MagickSpreadImage;
    np = 1;
  }
  else if (strstr(method,"strip"))
  {
    have_ch = 0;
    func    = &MagickStripImage;
    np = 0;
  }
  else if (strstr(method,"swirl"))
  {
    // parameters: radius
    have_ch = 0;
    func    = &MagickSwirlImage;
    np = 1;
  }
  else if (strstr(method,"thumbnail"))
  {
    // parameters: radius
    have_ch = 0;
    func    = &_MagickThumbnailImage;
    np = 2;
  }
  else if (strstr(method,"transparentpaint"))
  {
    have_ch = 0;
    func    = &MagickTransparentPaintImage;
    np = -1;
  }
  else if (strstr(method,"transpose"))
  {
    have_ch = 0;
    func    = &MagickTransposeImage;
    np = 0;
  }
  else if (strstr(method,"transverse"))
  {
    have_ch = 0;
    func    = &MagickTransverseImage;
    np = 0;
  }
  else if (strstr(method,"trim"))
  {
    // parameters: fuzz
    have_ch = 0;
    func    = &MagickTrimImage;
    np = 1;
  }
  else if (strstr(method,"uniquecolors"))
  {
    have_ch = 0;
    func    = &MagickUniqueImageColors;
    np = 0;
  }
  else if (strstr(method,"unsharpmask"))
  {
    // parameters: radius, sigma, amount, threshold
    have_ch = 1;
    func    = &MagickUnsharpMaskImage;
    func_ch = &MagickUnsharpMaskImageChannel;
    np = 4;
  }
  else if (strstr(method,"vignette"))
  {
    have_ch = 0;
    func = &_MagickVignetteImage;
    np = 4;
  }
  else if (strstr(method,"wave"))
  {
    have_ch = 0;
    func = &MagickWaveImage;
    np = 2;
  }
  else if (strstr(method,"whitethreshold"))
  {
    have_ch = 0;
    func = &MagickWhiteThresholdImage;
    np = -1;
  }
  else
  {
    fprintf(stderr, "Method %s not recognized!\n", method);
    print_recognized_methods();
    rerror(THIS_SOLVER  ": Cannot continue. Method is missing!\n");
  }

  //
  // process channels if necessary
  //
  if (have_ch)
  {
    lb = strstr(method,"[");
    rb = strstr(method,"]");
    if (lb && rb)
    {
      ar = (strstr(lb, "r") != 0);  // RedChannel
      ac = (strstr(lb, "c") != 0);  // CyanChannel
      ag = (strstr(lb, "g") != 0);  // GreenChannel
      am = (strstr(lb, "m") != 0);  // MagentaChannel
      ab = (strstr(lb, "b") != 0);  // BlueChannel
      ay = (strstr(lb, "y") != 0);  // YellowChannel
      ao = (strstr(lb, "o") != 0);  // OpacityChannel
      aB = (strstr(lb, "B") != 0);  // BlackChannel
      aM = (strstr(lb, "M") != 0);  // MatteChannel
    }
  }

  //
  // first process functions with special arguments
  // (that cannot be casted to double * for arguments)
  //
  if(func == &MagickAddNoiseImage)
  {
    if(nargs != np+1 && nargs != np)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_STRING)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires string as a second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    c = class_char_pointer(en);
    switch(*c)
    {
      case 'u':
      case 'U':
        noise_type = UniformNoise;
        break;

      case 'g':
      case 'G':
      case 'n':
      case 'N':
        noise_type = GaussianNoise;
        break;

      case 'm':
      case 'M':
        noise_type = MultiplicativeGaussianNoise;
        break;

      case 'i':
      case 'I':
        noise_type = ImpulseNoise;
        break;

      case 'l':
      case 'L':
        noise_type = LaplacianNoise;
        break;

      case 'p':
      case 'P':
        noise_type = PoissonNoise;
        break;

      default:
        noise_type = UniformNoise;
        break;
    }
    ent_Clean(en);

    if (nargs==3)
    {
      en = bltin_get_ent (args[2]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      sprintf(default_text, "%g", class_double(en));
      ent_Clean(en);

      MagickSetImageArtifact(magick_wand[default_wand_idx], "attenuate", default_text);
    }

    if (ar)
      func_ch(magick_wand[default_wand_idx], RedChannel, noise_type);
    if (ac)
      func_ch(magick_wand[default_wand_idx], CyanChannel,noise_type);
    if (ag)
      func_ch(magick_wand[default_wand_idx], GreenChannel,noise_type);
    if (am)
      func_ch(magick_wand[default_wand_idx], MagentaChannel,noise_type);
    if (ab)
      func_ch(magick_wand[default_wand_idx], BlueChannel,noise_type);
    if (ay)
      func_ch(magick_wand[default_wand_idx], YellowChannel,noise_type);
    if (ao)
      func_ch(magick_wand[default_wand_idx], OpacityChannel,noise_type);
    if (aB)
      func_ch(magick_wand[default_wand_idx], BlackChannel,noise_type);
    if (aM)
      func_ch(magick_wand[default_wand_idx], MatteChannel,noise_type);
    if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
      func(magick_wand[default_wand_idx], noise_type);
  }
  else if (func == &MagickResizeImage)
  {
    // special methods:
    if(nargs != 4)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x2 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires two entries vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=2)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires two entries vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = mdrV0(mdr_param,0);
    params[1] = mdrV0(mdr_param,1);
    ent_Clean(en);

    // third parameter is the filter name
    en = bltin_get_ent (args[2]);
    if (ent_type(en) != MATRIX_DENSE_STRING)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires string as a second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    FilterTypes filter;
    c = class_char_pointer(en);
    if (strstr(c, "bes"))
      filter = BesselFilter;
    else if (strstr(c, "bl"))
      filter = BlackmanFilter;
    else if (strstr(c,"box"))
      filter = BoxFilter;
    else if (strstr(c,"cat"))
      filter = CatromFilter;
    else if (strstr(c,"cub"))
      filter = CubicFilter;
    else if (strstr(c,"han"))
      filter = HanningFilter;
    else if (strstr(c,"her"))
      filter = HermiteFilter;
    else if (strstr(c,"lan"))
      filter = LanczosFilter;
    else if (strstr(c,"mit"))
      filter = MitchellFilter;
    else if (strstr(c,"quad"))
      filter = QuadraticFilter;
    else if (strstr(c,"sinc"))
      filter = SincFilter;
    else if (strstr(c,"tri"))
      filter = TriangleFilter;
    else
      filter = LanczosFilter;

    ent_Clean(en);

    // fourth is the blur
    en = bltin_get_ent (args[3]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires real 'blur' as fourth argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[2] = class_double(en);
    ent_Clean(en);

    func(magick_wand[default_wand_idx], (size_t) params[0], (size_t) params[1], filter, params[2]);
  }
  else if (func == &MagickRotateImage)
  {
    if(nargs != 3)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=3 && mdr_param->nrow * mdr_param->ncol !=4)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = mdrV0(mdr_param,0);
    params[1] = mdrV0(mdr_param,1);
    params[2] = mdrV0(mdr_param,2);
    if (mdr_param->nrow * mdr_param->ncol ==4 )
      params[3] = MdrV0(mdr_param,3);
    else
      params[3] = 1;
    ent_Clean(en);

    // fourth is the blur
    en = bltin_get_ent (args[2]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[4] = class_double(en);
    ent_Clean(en);

    PixelWand *pw = NewPixelWand();
    PixelSetRed(pw, params[0]);
    PixelSetGreen(pw, params[1]);
    PixelSetBlue(pw, params[2]);
    PixelSetOpacity(pw, params[3]);
    MagickRotateImage(magick_wand[default_wand_idx], pw,params[4]);
    pw = DestroyPixelWand(pw);
  }
  else if (func == &MagickShearImage)
  {
    if(nargs != 4)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=3 && mdr_param->nrow * mdr_param->ncol !=4)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = mdrV0(mdr_param,0);
    params[1] = mdrV0(mdr_param,1);
    params[2] = mdrV0(mdr_param,2);
    if (mdr_param->nrow * mdr_param->ncol ==4 )
      params[3] = MdrV0(mdr_param,3);
    else
      params[3] = 1;
    ent_Clean(en);

    // third is x_shear
    en = bltin_get_ent (args[2]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[4] = class_double(en);
    ent_Clean(en);

    // fourth is y_shear
    en = bltin_get_ent (args[3]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[5] = class_double(en);
    ent_Clean(en);

    PixelWand *pw = NewPixelWand();
    PixelSetRed(pw, params[0]);
    PixelSetGreen(pw, params[1]);
    PixelSetBlue(pw, params[2]);
    PixelSetOpacity(pw, params[3]);

    MagickShearImage(magick_wand[default_wand_idx], pw, params[4], params[5]);

    pw = DestroyPixelWand(pw);
  }
  else if (func == &MagickTransparentPaintImage)
  {
    if(nargs != 5)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
      en = bltin_get_ent (args[1]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      mdr_param = ent_data(en);
      if (mdr_param->nrow * mdr_param->ncol !=3 && mdr_param->nrow * mdr_param->ncol !=4)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[0] = mdrV0(mdr_param,0);
      params[1] = mdrV0(mdr_param,1);
      params[2] = mdrV0(mdr_param,2);
      if (mdr_param->nrow * mdr_param->ncol ==4 )
        params[3] = MdrV0(mdr_param,3);
      else
        params[3] = 1;
      ent_Clean(en);

      // third is alpha
      en = bltin_get_ent (args[2]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[4] = class_double(en);
      ent_Clean(en);

      // fourth is fuzz
      en = bltin_get_ent (args[3]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[5] = class_double(en);
      ent_Clean(en);

      // fifth is invert
      en = bltin_get_ent (args[4]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[6] = class_double(en);
      ent_Clean(en);

      PixelWand *pw = NewPixelWand();
      PixelSetRed(pw, params[0]);
      PixelSetGreen(pw, params[1]);
      PixelSetBlue(pw, params[2]);
      PixelSetOpacity(pw, params[3]);

      MagickTransparentPaintImage(magick_wand[default_wand_idx], pw, params[4], params[5],
                                  (MagickBooleanType) params[6]);

      pw = DestroyPixelWand(pw);
  }
  else if (func == &MagickWhiteThresholdImage)
  {
    if(nargs != 2)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=3 && mdr_param->nrow * mdr_param->ncol !=4)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = mdrV0(mdr_param,0);
    params[1] = mdrV0(mdr_param,1);
    params[2] = mdrV0(mdr_param,2);
    if (mdr_param->nrow * mdr_param->ncol ==4 )
      params[3] = MdrV0(mdr_param,3);
    else
      params[3] = 1;
    ent_Clean(en);

    PixelWand *pw = NewPixelWand();
    PixelSetRed(pw, params[0]);
    PixelSetGreen(pw, params[1]);
    PixelSetBlue(pw, params[2]);
    PixelSetOpacity(pw, params[3]);

    MagickWhiteThresholdImage(magick_wand[default_wand_idx], pw);

    pw = DestroyPixelWand(pw);
  }
  else if (func == &MagickBlackThresholdImage)
  {
    if(nargs != 2)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=3 && mdr_param->nrow * mdr_param->ncol !=4)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = mdrV0(mdr_param,0);
    params[1] = mdrV0(mdr_param,1);
    params[2] = mdrV0(mdr_param,2);
    if (mdr_param->nrow * mdr_param->ncol ==4 )
      params[3] = mdrV0(mdr_param,3);
    else
      params[3] = 1;
    ent_Clean(en);

    PixelWand *pw = NewPixelWand();
    PixelSetRedQuantum(pw, (Quantum) params[0]);
    PixelSetGreenQuantum(pw, (Quantum) params[1]);
    PixelSetBlueQuantum(pw, (Quantum) params[2]);
    PixelSetOpacityQuantum(pw, (Quantum) params[3]);

    MagickBlackThresholdImage(magick_wand[default_wand_idx], pw);

    pw = DestroyPixelWand(pw);
  }
  else if (func == &MagickAnnotateImage)
  {
    if(nargs != 4)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // first argument: text
    char *ann=0;
    Ent *es = bltin_get_ent (args[1]);
    if (ent_type(es) != MATRIX_DENSE_STRING)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires text as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    ann = class_char_pointer(es);

    // second argument: [x,y]
    en = bltin_get_ent (args[2]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires [x,y] real vector as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=2)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires [x,y] real vector as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = mdrV0(mdr_param,0);
    params[1] = mdrV0(mdr_param,1);
    ent_Clean(en);

    // third argument: angle
    en = bltin_get_ent (args[3]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires real scalar 'angle' as fourth argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[2] = class_double(en);
    ent_Clean(en);

    DrawingWand *drawing_wand = NewDrawingWand();
    rval = MagickAnnotateImage(magick_wand[default_wand_idx],drawing_wand,
                               params[0],params[1],params[2],ann);
    drawing_wand = DestroyDrawingWand(drawing_wand);
    ent_Clean(es);
  }
  else if (func == &MagickColorDecisionListImage)
  {
    char *filename=0;

    if(nargs != 2 && nargs!=5)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    if (nargs == 2)
    {
      Ent *en = bltin_get_ent (args[1]);
      if (ent_type(en) != MATRIX_DENSE_STRING)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires filename as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      filename = class_char_pointer(en);
      if (!filename)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires filename as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      FILE *fr = fopen (filename, "r");
      if (!fr)
      {
        fprintf(stderr, THIS_SOLVER  ": file %s does not exist!\n", filename);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      i=0;
      while( ( default_text[i] = fgetc(fr) ) != EOF  && i<MAX_LEN_DEFAULT_TEXT)
        i++;
      if (i==MAX_LEN_DEFAULT_TEXT)
      {
        fprintf(stderr, THIS_SOLVER  ": file %s is too big! It should contain a single profile!\n", filename);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      default_text[i]='\0';
      fclose(fr);
      ent_Clean(en);
    }
    else
    {
      // we create XML string from the user provided input data

      // second argument is a 1x3 matrix
      en = bltin_get_ent (args[1]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'offset' real vector as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      mdr_param = ent_data(en);
      if (mdr_param->nrow * mdr_param->ncol !=3)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'offset' real vector as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[0] = mdrV0(mdr_param,0);
      params[1] = mdrV0(mdr_param,1);
      params[2] = mdrV0(mdr_param,2);
      ent_Clean(en);

      // third argument is a 1x3 matrix
      en = bltin_get_ent (args[2]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'slope' real vector as third argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      mdr_param = ent_data(en);
      if (mdr_param->nrow * mdr_param->ncol !=3)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'slope' real vector as third argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[3] = mdrV0(mdr_param,0);
      params[4] = mdrV0(mdr_param,1);
      params[5] = mdrV0(mdr_param,2);
      ent_Clean(en);

      // fourth argument is a 1x3 matrix
      en = bltin_get_ent (args[3]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'power' real vector as fourth argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      mdr_param = ent_data(en);
      if (mdr_param->nrow * mdr_param->ncol !=3)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'power' real vector as fourth argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[6] = mdrV0(mdr_param,0);
      params[7] = mdrV0(mdr_param,1);
      params[8] = mdrV0(mdr_param,2);
      ent_Clean(en);

      // fifth argument is a scalar
      en = bltin_get_ent (args[4]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'saturation' real scalar as fifth argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[9] = class_double(en);
      ent_Clean(en);

      // now assemble XML file
      i=0;
      char *current_text=default_text;
      //
      sprintf(&current_text[i],"<ColorCorrectionCollection xmlns=\"urn:ASC:CDL:v0.01\">\n");
      i = strlen(current_text);
      //
      sprintf(&current_text[i],"<ColorCorrection id=\"rlab0000\">\n");
      i = strlen(current_text);
      //
      sprintf(&current_text[i],"<SOPNode>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"<Description>RLaB was here!</Description>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"<Offset> %g %g %g </Offset>\n",params[0],params[1],params[2]);
      i = strlen(current_text);
      sprintf(&current_text[i],"<Slope> %g %g %g </Slope>\n",params[3],params[4],params[5]);
      i = strlen(current_text);
      sprintf(&current_text[i],"<Power> %g %g %g </Power>\n",params[6],params[7],params[8]);
      i = strlen(current_text);
      sprintf(&current_text[i],"</SOPNode>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"<SatNode>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"<Saturation> %g </Saturation>\n",params[9]);
      i = strlen(current_text);
      sprintf(&current_text[i],"</SatNode>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"</ColorCorrection>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"</ColorCorrectionCollection>\n");
      i = strlen(current_text);
      current_text[i] = '\0';
    }

    if (im_debug)
      fprintf(stderr, "%s", default_text);

    bval = MagickColorDecisionListImage(magick_wand[default_wand_idx], default_text);
    if (bval == MagickFalse)
      WriteWandException(magick_wand[default_wand_idx]);
  }
  else if (func == &MagickColorizeImage)
  {
    if(nargs != 2)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 or 1x4 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'color' real 3- or 4-vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow *mdr_param->ncol != 3 && mdr_param->nrow *mdr_param->ncol != 4)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'color' real 3- or 4-vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }

    PixelWand *c = NewPixelWand();
    PixelSetRed   (c, MdrV0(mdr_param,0));
    PixelSetGreen (c, MdrV0(mdr_param,1));
    PixelSetBlue  (c, MdrV0(mdr_param,2));

    PixelWand *o = NewPixelWand();
    if (mdr_param->nrow *mdr_param->ncol == 4)
      PixelSetOpacity(o, MdrV0(mdr_param,3));
    else
      PixelSetOpacity(o, 1.0);

    MagickColorizeImage(magick_wand[default_wand_idx], c, o);

    c = DestroyPixelWand(c);
    o = DestroyPixelWand(o);
  }
  else if (func == &MagickColorMatrixImage)
  {
    if(nargs != 2)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow != mdr_param->ncol)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    if (mdr_param->nrow != 3 && mdr_param->nrow != 4 && mdr_param->nrow != 5)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }

    // construct string: rows in the matrix
    i=0;
    char *current_text=default_text;

    for (i=0; i<mdr_param->nrow; i++)
      for (j=0; j<mdr_param->ncol; j++)
      {
        sprintf(&current_text[i],"%g ", Mdr0(mdr_param,i,j));
        i = strlen(current_text);
      }
    current_text[i-1] ='\0';  // remove extraspace at the end

    KernelInfo *mykernel = AcquireKernelInfo(default_text);
    MagickColorMatrixImage(magick_wand[default_wand_idx],mykernel);

    mykernel = DestroyKernelInfo(mykernel);
    ent_Clean(en);
  }
  else if  (func == &MagickConvolveImage)
  {
    if(nargs != 2)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    MDR *order=0;

    // second argument is a square matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow != mdr_param->ncol)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    order = mdr_Transpose(mdr_param);

    if (ar)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], RedChannel, (size_t) MNR(order), MDRPTR(order));
    if (ac)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], CyanChannel, (size_t) MNR(order), MDRPTR(order));
    if (ag)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], GreenChannel, (size_t) MNR(order), MDRPTR(order));
    if (am)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], MagentaChannel, (size_t) MNR(order), MDRPTR(order));
    if (ab)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], BlueChannel, (size_t) MNR(order), MDRPTR(order));
    if (ay)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], YellowChannel, (size_t) MNR(order), MDRPTR(order));
    if (ao)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], OpacityChannel, (size_t) MNR(order), MDRPTR(order));
    if (aB)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], BlackChannel, (size_t) MNR(order), MDRPTR(order));
    if (aM)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], MatteChannel, (size_t) MNR(order), MDRPTR(order));
    if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
      MagickConvolveImage(magick_wand[default_wand_idx], (size_t) MNR(order), MDRPTR(order));

    mdr_Destroy(order);
    ent_Clean(en);
  }
  else
  {
    //
    // generic method with generic parameters
    //
    if (nargs==2 && np>1)
    {
      // all parameters are in the single array
      en = bltin_get_ent (args[1]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

      mdr_param = ent_data(en);
      if (mdr_param->nrow * mdr_param->ncol!=np)
        rerror(THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      for (i=0; i<np; i++)
        params[i] = mdrV0(mdr_param, i);

      ent_Clean(en);
    }
    else if (nargs==(np+1) && np>1)
    {
      // read parameters from arguments to the function
      for (i=1; i<=np; i++)
      {
        en = bltin_get_ent (args[i]);
        if (ent_type(en) == MATRIX_DENSE_REAL)
          params[i-1] = class_double(en);
        ent_Clean(en);
      }
    }
    else if ((np+1)!=nargs)
      rerror(THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    //
    // generic method with generic parameters
    //
    switch (np)
    {
      case 0:
        // no parameter functions
        if (!have_ch)
        {
          func(magick_wand[default_wand_idx]);
        }
        else
        {
          if (ar)
            func_ch(magick_wand[default_wand_idx], RedChannel);
          if (ac)
            func_ch(magick_wand[default_wand_idx], CyanChannel);
          if (ag)
            func_ch(magick_wand[default_wand_idx], GreenChannel);
          if (am)
            func_ch(magick_wand[default_wand_idx], MagentaChannel);
          if (ab)
            func_ch(magick_wand[default_wand_idx], BlueChannel);
          if (ay)
            func_ch(magick_wand[default_wand_idx], YellowChannel);
          if (ao)
            func_ch(magick_wand[default_wand_idx], OpacityChannel);
          if (aB)
            func_ch(magick_wand[default_wand_idx], BlackChannel);
          if (aM)
            func_ch(magick_wand[default_wand_idx], MatteChannel);
          if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
            func(magick_wand[default_wand_idx]);
        }
        break;

      case 1:
        if (!have_ch)
        {
          func(magick_wand[default_wand_idx], params[0]);
        }
        else
        {
          if (ar)
            func_ch(magick_wand[default_wand_idx], RedChannel, params[0]);
          if (ac)
            func_ch(magick_wand[default_wand_idx], CyanChannel,params[0]);
          if (ag)
            func_ch(magick_wand[default_wand_idx], GreenChannel,params[0]);
          if (am)
            func_ch(magick_wand[default_wand_idx], MagentaChannel,params[0]);
          if (ab)
            func_ch(magick_wand[default_wand_idx], BlueChannel,params[0]);
          if (ay)
            func_ch(magick_wand[default_wand_idx], YellowChannel,params[0]);
          if (ao)
            func_ch(magick_wand[default_wand_idx], OpacityChannel,params[0]);
          if (aB)
            func_ch(magick_wand[default_wand_idx], BlackChannel,params[0]);
          if (aM)
            func_ch(magick_wand[default_wand_idx], MatteChannel,params[0]);
          if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
            func(magick_wand[default_wand_idx], params[0]);
        }
        break;

      case 2:
        // generic two parameter functions
        if (!have_ch)
        {
          func(magick_wand[default_wand_idx], params[0], params[1]);
        }
        else
        {
          if (ar)
            func_ch(magick_wand[default_wand_idx], RedChannel,params[0], params[1]);
          if (ac)
            func_ch(magick_wand[default_wand_idx], CyanChannel,params[0], params[1]);
          if (ag)
            func_ch(magick_wand[default_wand_idx], GreenChannel,params[0], params[1]);
          if (am)
            func_ch(magick_wand[default_wand_idx], MagentaChannel,params[0], params[1]);
          if (ab)
            func_ch(magick_wand[default_wand_idx], BlueChannel,params[0], params[1]);
          if (ay)
            func_ch(magick_wand[default_wand_idx], YellowChannel,params[0], params[1]);
          if (ao)
            func_ch(magick_wand[default_wand_idx], OpacityChannel,params[0], params[1]);
          if (aB)
            func_ch(magick_wand[default_wand_idx], BlackChannel,params[0], params[1]);
          if (aM)
            func_ch(magick_wand[default_wand_idx], MatteChannel,params[0], params[1]);
          if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
            func(magick_wand[default_wand_idx], params[0], params[1]);
        }
        break;

      case 3:
        // three parameter functions
        if (!have_ch)
          func(magick_wand[default_wand_idx], params[0], params[1],params[2]);
        else
        {
          if (ar)
            func_ch(magick_wand[default_wand_idx],RedChannel,params[0],params[1],params[2]);
          if (ac)
            func_ch(magick_wand[default_wand_idx],CyanChannel,params[0],params[1],params[2]);
          if (ag)
            func_ch(magick_wand[default_wand_idx],GreenChannel,params[0],params[1],params[2]);
          if (am)
            func_ch(magick_wand[default_wand_idx],MagentaChannel,params[0],params[1],params[2]);
          if (ab)
            func_ch(magick_wand[default_wand_idx],BlueChannel,params[0],params[1],params[2]);
          if (ay)
            func_ch(magick_wand[default_wand_idx],YellowChannel,params[0],params[1],params[2]);
          if (ao)
            func_ch(magick_wand[default_wand_idx],OpacityChannel,params[0],params[1],params[2]);
          if (aB)
            func_ch(magick_wand[default_wand_idx],BlackChannel,params[0],params[1],params[2]);
          if (aM)
            func_ch(magick_wand[default_wand_idx],MatteChannel,params[0],params[1],params[2]);
          if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
            func(magick_wand[default_wand_idx], params[0], params[1],params[2]);
        }
        break;

      case 4:
        if (func == &MagickChopImage || func == &MagickCropImage || func == &MagickExtentImage)
        {
          func(magick_wand[default_wand_idx], (size_t) params[0], (size_t) params[1],
               (size_t) params[2], (size_t) params[3]);
        }
        else if (func == &MagickLiquidRescaleImage)
        {
          func(magick_wand[default_wand_idx], (size_t) params[0], (size_t) params[1], params[2],params[3]);
        }
        break;

      case 5:
        break;

      case 6:
        if (func == &MagickAffineTransformImage)
        {
          fprintf(stderr, "'affinetransform' is not working!\n");
          AffineMatrix aff_mat;
          DrawingWand *aff_dw=NewDrawingWand();
          aff_mat.sx = params[0];
          aff_mat.rx = params[1];
          aff_mat.ry = params[2];
          aff_mat.sy = params[3];
          aff_mat.tx = params[4];
          aff_mat.ty = params[5];
          DrawAffine(aff_dw, &aff_mat);
//           DrawRender(aff_dw);
          if (func(magick_wand[default_wand_idx], aff_dw) == MagickFalse)
            WriteWandException(magick_wand[default_wand_idx]);
//           DrawRender(aff_dw);
//           MagickDrawImage(magick_wand[default_wand_idx],aff_dw);
          aff_dw = DestroyDrawingWand(aff_dw);
        }
        else if (func == &MagickBorderImage)
        {
          PixelWand *pw = NewPixelWand();
          PixelSetRed     (pw, params[0]);
          PixelSetGreen   (pw, params[1]);
          PixelSetBlue    (pw, params[2]);
          PixelSetOpacity (pw, params[3]);
          func(magick_wand[default_wand_idx], pw, (size_t) params[4], (size_t) params[5]);
          pw = DestroyPixelWand(pw);
        }
        break;

      default:
        break;
    }
  }

  ent_Clean(e1);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (rval);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}


//
// set default wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".wand"
Ent *
ent_im_wand (int nargs, Datum args[])
{
  Ent *e1=0, *rent, *R;
  int i,k,icurr;
  short int *image_idx = (short int *) default_text;
  Btree *bw=0;

  if (nargs != 1 && nargs != 0)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Get/set default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "([wand]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'wand' is an index of the wand.\n");
    rerror ( THIS_SOLVER ": requires 1 argument");
  }

  // special case:
  //  no arguments and library not initialized
  //  just return empty list
  if (nargs ==0 && !InitMagickWandLib)
  {
    rent = ent_Create ();
    ent_data (rent) = btree_Create();;
    ent_type (rent) = BTREE;
    return rent;
  }


  if(!InitMagickWandLib)
  {
    MagickWandGenesis();
    InitMagickWandLib = 1;
  }

  if (nargs==1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) != MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": First argument 'wand' must be an integer!");
    i = (int) class_double(e1) - 1;
    if (i>=0 && i<MAX_NUMBER_MAGICK_WANDS)
      default_wand_idx = i;
  }
  if (!magick_wand[default_wand_idx])
    magick_wand[default_wand_idx] = NewMagickWand();

  ent_Clean(e1);

  bw = btree_Create();

  // wand:
  R = ent_Create();
  ent_data (R) = mdi_CreateScalar(default_wand_idx + 1);
  ent_type (R) = MATRIX_DENSE_REAL;
  install  (bw, "wand", R);

  // figure out the active image
  icurr = MagickGetIteratorIndex(magick_wand[default_wand_idx]);
  if (icurr+1)
    install  (bw, "image", ent_Create_Rlab_Int(icurr));
  else
    install  (bw, "image", ent_Assign_Rlab_MDR(NULL));

  // figure out all the images in the wand
  MDR *imgs=0;
  MagickResetIterator(magick_wand[default_wand_idx]);
  k = -1;
  while(MagickNextImage(magick_wand[default_wand_idx]))
    image_idx[++k] = MagickGetIteratorIndex(magick_wand[default_wand_idx]);
  k++;
  if (k)
  {
    imgs = mdi_Create(1,k);
    for (i=0;i<k;i++)
      MdiV0(imgs,i) = image_idx[i] + 1;
    // put back the original image on the wand
    MagickSetIteratorIndex(magick_wand[default_wand_idx],icurr);
  }
  else
    imgs = mdi_Create(0,0);

  install  (bw, "images", ent_Assign_Rlab_MDR(imgs));

  return ent_Assign_Rlab_BTREE(bw);
}


#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".write"
Ent *
ent_im_wand_write_image (int nargs, Datum args[])
{
  Ent *e1=0, *rent;
  char *filename=0;

  int rval,save_all=0;

  if (nargs > 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Write active image or all images from default Magick Wand to file.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "([filename],[all]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'filename' is where active image will be saved, while\n");
    fprintf (stdout,
             THIS_SOLVER ":   'all' is an indicator to save all images, and not just active image.\n");
    fprintf (stdout,
             THIS_SOLVER ": See manual for details.\n");
    rerror ( THIS_SOLVER ": requires 2 or 3 arguments");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) != MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": Second argument 'filename' must be string!\n");
    filename  = class_char_pointer(e1);
    if (!filename)
      rerror(THIS_SOLVER ": Second argument 'filename' must be string!\n");
  }

  if (nargs == 2)
    save_all = 1;

  if (save_all)
    rval = MagickWriteImages(magick_wand[default_wand_idx], filename, MagickTrue);
  else
    rval = MagickWriteImage (magick_wand[default_wand_idx], filename);
  if (rval == MagickFalse)
  {
    WriteWandException(magick_wand[default_wand_idx]);
    rerror(THIS_SOLVER ": Operation failed!");
  }

  ent_Clean(e1);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar(0);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

//
// get/set wand artifacts
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".art"
Ent *
ent_im_artifacts (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  int j;
  size_t n=100;

  char *pattern_all = "*";
  char *pattern_mine=0;
  char **c, *opt=0, *new_val=0;
  double dummy;

  MDS *s=0;

  if (nargs > 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": List/set/get artifacts of default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "([artifacts],[value]),\n");
    rerror ( THIS_SOLVER ": requires none, 1 or 2 arguments");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs>=1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_STRING)
      pattern_mine  = class_char_pointer(e1);
  }

  if (nargs==2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) == MATRIX_DENSE_STRING)
      new_val = class_char_pointer(e2);
    else if (ent_type(e2) == MATRIX_DENSE_REAL)
    {
      dummy = class_double(e2);
      sprintf(default_text,"%g", dummy);
      new_val = default_text;
    }
  }

  if (nargs==0 || !pattern_mine)
  {
    c =  (char **) MagickGetImageArtifacts(magick_wand[default_wand_idx], pattern_all, &n);
    s = mds_Create(1*(n>0), n);
    if (n)
    {
      for (j=0; j<n; j++)
      {
        MdsV0(s,j) = cpstr(c[j]);
      }
    }
    MagickRelinquishMemory(c);
  }
  else if (nargs==1 && pattern_mine)
  {
    opt = MagickGetImageArtifact(magick_wand[default_wand_idx], pattern_mine);
    if (opt)
      s = mds_CreateScalar(opt);
    else
      s = mds_Create(0, 0);
    MagickRelinquishMemory(opt);
  }
  else if (nargs==2 && pattern_mine && new_val)
  {
    MagickSetImageArtifact(magick_wand[default_wand_idx], pattern_mine, new_val);
    opt = MagickGetImageArtifact(magick_wand[default_wand_idx], pattern_mine);
    if (opt)
      s = mds_CreateScalar( opt );
    else
      s = mds_Create(0, 0);
    MagickRelinquishMemory(opt);
  }
  else if (nargs==2 && pattern_mine && ent_type(e2)==UNDEF)
  {
    opt = MagickGetImageArtifact(magick_wand[default_wand_idx], pattern_mine);
    if (opt)
      MagickDeleteImageArtifact(magick_wand[default_wand_idx], pattern_mine);
    s = mds_Create(0, 0);
    MagickRelinquishMemory(opt);
  }
  else
    fprintf(stderr, THIS_SOLVER  ": Help! I don't know what I need to do!\n");

  ent_Clean(e1);
  ent_Clean(e2);

  rent = ent_Create ();
  ent_data (rent) = s;
  ent_type (rent) = MATRIX_DENSE_STRING;
  return (rent);
}


//
// get/set wand options
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".opt"
Ent *
ent_im_options (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  int j;
  size_t n=100;

  char *pattern_all = "*";
  char *pattern_mine=0;
  char **c, *opt=0, *new_val=0;
  double dummy;

  MDS *s=0;

  if (nargs > 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": List/set/get options of default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "([option],[value]),\n");
    rerror ( THIS_SOLVER ": requires none, 1 or 2 arguments");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs>=1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_STRING)
      pattern_mine  = class_char_pointer(e1);
  }

  if (nargs==2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) == MATRIX_DENSE_STRING)
      new_val = class_char_pointer(e2);
    else if (ent_type(e2) == MATRIX_DENSE_REAL)
    {
      dummy = class_double(e2);
      sprintf(default_text,"%g", dummy);
      new_val = default_text;
    }
  }

  if (nargs==0 || !pattern_mine)
  {
    c =  (char **) MagickGetOptions(magick_wand[default_wand_idx], pattern_all, &n);
    s = mds_Create(1*(n>0), n);
    if (n)
    {
      for (j=0; j<n; j++)
      {
        MdsV0(s,j) = cpstr(c[j]);
      }
    }
    MagickRelinquishMemory(c);
  }
  else if (nargs==1 && pattern_mine)
  {
    opt = MagickGetOption(magick_wand[default_wand_idx], pattern_mine);
    if (opt)
      s = mds_CreateScalar( opt );
    else
      s = mds_Create(0, 0);
    MagickRelinquishMemory(opt);
  }
  else if (nargs==2 && pattern_mine && new_val)
  {
    MagickSetOption(magick_wand[default_wand_idx], pattern_mine, new_val);
    opt = MagickGetOption(magick_wand[default_wand_idx], pattern_mine);
    if (opt)
      s = mds_CreateScalar( opt );
    else
      s = mds_Create(0, 0);
    MagickRelinquishMemory(opt);
  }
  else if (nargs==2 && pattern_mine && ent_type(e2)==UNDEF)
  {
    opt = MagickGetOption(magick_wand[default_wand_idx], pattern_mine);
    if (opt)
      MagickDeleteOption(magick_wand[default_wand_idx], pattern_mine);
    s = mds_Create(0, 0);
    MagickRelinquishMemory(opt);
  }
  else
    fprintf(stderr, THIS_SOLVER  ": Help! I don't know what I need to do!\n");

  ent_Clean(e1);
  ent_Clean(e2);

  rent = ent_Create ();
  ent_data (rent) = s;
  ent_type (rent) = MATRIX_DENSE_STRING;
  return (rent);
}


//
// display an image
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".disp"
#include <unistd.h>
Ent *
ent_im_display_wand (int nargs, Datum args[])
{
  Ent *e1=0;
  char *xsrv;
  char *c=0;

  if (nargs > 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Display active image from default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "([xserv]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'xsrv' is the name of x-server.\n");
    fprintf (stdout,
             THIS_SOLVER ": To return the control to RLAB, the display must be closed. See manual for details.\n");
    rerror ( THIS_SOLVER ": requires none or 1 argument");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_STRING)
      c  = class_char_pointer(e1);
  }

  // figure out the X-server
  if (c)
    xsrv = c;
  else
    xsrv = getenv("DISPLAY");

  pid_t pid = vfork();
  if(pid == 0)
  {
    MagickDisplayImage(magick_wand[default_wand_idx], xsrv);
    exit(0);              // child exits
  }

//   if (rval == MagickFalse)
//   {
//     WriteWandException(magick_wand[default_wand_idx]);
//     rerror(THIS_SOLVER ": Operation failed!");
//   }

  ent_Clean(e1);
  return ent_Create_Rlab_Success();
}


//
// load an image into the wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".read"
Ent *
ent_im_read_image_to_wand (int nargs, Datum args[])
{
  Ent *e1=0,  *e2=0;
  char *filename=0, *c=0;
  ssize_t i;

  if (nargs != 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Read image from file into default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(filename[,pos]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'filename', contains the image(s), while\n");
    fprintf (stdout,
             THIS_SOLVER ":   'pos', contains the position in wand forthwith from image are loaded.\n");
    fprintf (stdout,
             THIS_SOLVER ": See manual for details.\n");
    rerror ( THIS_SOLVER ": requires 1 or 2 arguments");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror(THIS_SOLVER ": First argument 'filename' must be string!");
  filename  = class_char_pointer(e1);
  if (!filename)
    rerror(THIS_SOLVER ": First argument 'filename' must be string!");

  if (nargs == 1)
  {
    MagickSetLastIterator(magick_wand[default_wand_idx]);
  }
  else
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) == MATRIX_DENSE_STRING)
    {
      c  = class_char_pointer(e1);
      if (c)
      {
        switch (*c)
        {
          case 'c':
          case 'C':
            i = MagickGetIteratorIndex(magick_wand[default_wand_idx]);
            MagickSetIteratorIndex(magick_wand[default_wand_idx], i);
            break;

          case 'f':
          case 'F':
            MagickSetFirstIterator(magick_wand[default_wand_idx]);
            break;

          default:
            // case 'l':
            MagickSetLastIterator(magick_wand[default_wand_idx]);
            break;
        }
      }
    }
  }

  if (MagickReadImage(magick_wand[default_wand_idx], filename) == MagickFalse)
  {
    WriteWandException(magick_wand[default_wand_idx]);
    rerror(THIS_SOLVER ": Operation failed!");
  }

  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Create_Rlab_Success();
}


//
// initialize the library
//
Ent *
ent_im_init (int nargs, Datum args[])
{
  if (!InitMagickWandLib)
  {
    MagickWandGenesis();
    InitMagickWandLib = 1;
  }

  return ent_Create_Rlab_Success();
}


//
// iterate through images in a wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".iter"
Ent *
ent_im_wand_iterate_images (int nargs, Datum args[])
{
  Ent *e1=0;
  int j, rval=0;
  char *op=0;

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs != 0 && nargs != 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Iterate through images in default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "([op]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'op' = 'first' 'last' 'next' 'prev' or NNN determine the operation\n");
    rerror ( THIS_SOLVER ": requires none or 1 argument");
  }


  if (nargs == 0)
  {
    rval = (int) MagickGetIteratorIndex(magick_wand[default_wand_idx]) + 1;
  }
  else if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_STRING)
    {
      op = class_char_pointer(e1);
      switch(*op)
      {
        case 'f':
        case 'F':
          MagickSetFirstIterator(magick_wand[default_wand_idx]);
          rval=0;
          break;

        case 'n':
        case 'N':
          if (MagickNextImage(magick_wand[default_wand_idx]) == MagickFalse)
          {
            MagickPreviousImage(magick_wand[default_wand_idx]);
            rval=1;
          }
          else
            rval=0;
          break;

        case 'l':
        case 'L':
          MagickSetLastIterator(magick_wand[default_wand_idx]);
          rval=0;
          break;

        case 'p':
        case 'P':
          if (MagickPreviousImage(magick_wand[default_wand_idx]) == MagickFalse)
          {
            MagickNextImage(magick_wand[default_wand_idx]);
            rval=1;
          }
          else
            rval=0;
          break;

        default:
          rval = 2;
          break;
      }
    }
    else if (ent_type(e1) == MATRIX_DENSE_REAL)
    {
      j = (int) class_double(e1) - 1;
      if (MagickSetIteratorIndex(magick_wand[default_wand_idx], j) == MagickFalse)
        rval = 1;
      else
        rval = 0;
    }
    else
      rval = 4;
  }

  ent_Clean(e1);

  return ent_Create_Rlab_Double(rval);
}


//
// add image from one wand to another
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".add"
Ent *
ent_im_add (int nargs, Datum args[])
{
  Ent *e1=0;
  int j;

  if (nargs != 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Add images from Magick Wand to the default Magic Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(add_wand),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'add_wand' is the wand which images are added to the default one.\n");
    rerror ( THIS_SOLVER ": requires 1 argument");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": First argument 'add_wand' must be an integer!");
  j = (int) class_double(e1) - 1;
  if (!magick_wand[j])
    rerror(THIS_SOLVER ": First argument 'add_wand' must be an existing Magick Wand!");

  if (j>=0 && j<MAX_NUMBER_MAGICK_WANDS)
  {
    if(MagickAddImage(magick_wand[default_wand_idx], magick_wand[j])==MagickFalse)
    {
      WriteWandException(magick_wand[default_wand_idx]);
    }
  }

  ent_Clean(e1);

  return ent_Create_Rlab_Double(default_wand_idx+1);
}


//
// clone a wand and new wand becomes a default one
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".clone"
Ent *
ent_im_clone (int nargs, Datum args[])
{
  Ent *e1=0;
  int j;

  if (nargs != 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Clone default Magick Wand, and make new one the default.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(new_wand),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'new_wand' is clone of the default one, and the new default wand.\n");
    rerror ( THIS_SOLVER ": requires 1 argument");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": First argument 'wand' must be an integer!");
  j = (int) class_double(e1) - 1;

  if ((j>=0 && j<MAX_NUMBER_MAGICK_WANDS) && (default_wand_idx!=j))
  {
    // destroy or clear the target
    if (magick_wand[j])
      DestroyMagickWand(magick_wand[j]);

    // clone the first wand
    magick_wand[j] = CloneMagickWand(magick_wand[default_wand_idx]);
    default_wand_idx = j;
  }
  else
    fprintf(stderr,THIS_SOLVER  ": Index out of bounds!\n");

  ent_Clean(e1);
  return ent_Create_Rlab_Double(default_wand_idx+1);
}


//
// get/set wand properties
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".prop"
Ent *
ent_im_property (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  int j;
  size_t n=-1;

  char *pattern_all = "*";
  char *pattern_mine=0;
  char **c, *opt=0, *new_val=0;

  MDS *s=0;

  if (nargs > 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": List/set/get properties of a Magick Wand from all images.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "([property],value/]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'property' is the image property to be retrieved or set, while\n");
    fprintf (stdout,
             THIS_SOLVER ":   'value' is its value to be set.\n");
    rerror ( THIS_SOLVER ": requires 1 argument");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs>=1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_STRING)
      pattern_mine  = class_char_pointer(e1);
  }

  if (nargs==2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) == MATRIX_DENSE_STRING)
      new_val = class_char_pointer(e2);
  }

  if (nargs==0 || !pattern_mine)
  {
    c =  (char **) MagickGetImageProperties(magick_wand[default_wand_idx], pattern_all, &n);
    if (c)
    {
      if (n>0)
      {
        s = mds_Create(1, n);
        for (j=0; j<n; j++)
        {
          MdsV0(s,j) = cpstr(c[j]);
        }
      }
      MagickRelinquishMemory(c);
    }
  }
  else if (nargs==1 && pattern_mine)
  {
    opt = MagickGetImageProperty(magick_wand[default_wand_idx], pattern_mine);
    if (opt)
      s = mds_CreateScalar( opt );
    else
      s = mds_Create(0, 0);
    MagickRelinquishMemory(opt);
  }
  else if (nargs==2 && pattern_mine && new_val)
  {
    MagickSetImageProperty(magick_wand[default_wand_idx], pattern_mine, new_val);
    opt = MagickGetImageProperty(magick_wand[default_wand_idx], pattern_mine);
    if (opt)
      s = mds_CreateScalar( opt );
    else
      s = mds_Create(0, 0);
    MagickRelinquishMemory(opt);
  }
  else if (nargs==2 && pattern_mine && ent_type(e2)==UNDEF)
  {
    opt = MagickGetImageProperty(magick_wand[default_wand_idx], pattern_mine);
    if (opt)
      MagickDeleteImageProperty(magick_wand[default_wand_idx], pattern_mine);
    s = mds_Create(0, 0);
    MagickRelinquishMemory(opt);
  }
  else
  {
    s = mds_Create(0, 0);
    fprintf(stderr, THIS_SOLVER  ": Help! I don't know what I need to do!\n");
  }

  ent_Clean(e1);
  ent_Clean(e2);

  rent = ent_Create ();
  ent_data (rent) = s;
  ent_type (rent) = MATRIX_DENSE_STRING;
  return (rent);
}


//
// clear a wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".clear"
Ent *
ent_im_clear (int nargs, Datum args[])
{
  if (nargs)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Clear default Magick Wand from all images.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(),\n");
    rerror ( THIS_SOLVER ": requires no arguments");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  ClearMagickWand(magick_wand[default_wand_idx]);

  return ent_Create_Rlab_Int(default_wand_idx+1);
}


//
// close the default wand, if none is left close the library
//
Ent *
ent_im_close (int nargs, Datum args[])
{
  Ent *rent;
  int i;

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  magick_wand[default_wand_idx] = DestroyMagickWand(magick_wand[default_wand_idx]);
  magick_wand[default_wand_idx] = 0;
  default_wand_idx = -1;

  for(i=0; i<MAX_NUMBER_MAGICK_WANDS; i++)
  {
    if (magick_wand[i])
    {
      default_wand_idx = i;
      break;
    }
  }

  if (default_wand_idx == -1)
  {
    MagickWandTerminus();
    InitMagickWandLib = 0;
  }

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}


//
// exit the library close all wands
//
Ent *
ent_im_exit (int nargs, Datum args[])
{
  Ent *rent;
  int i;
  int rval=0;

  if (!InitMagickWandLib)
  {
    rval = 1;
    goto exit_exit;
  }

  default_wand_idx = -1;

  for(i=0; i<MAX_NUMBER_MAGICK_WANDS; i++)
  {
    if (magick_wand[i])
    {
      magick_wand[i] = DestroyMagickWand(magick_wand[i]);
      magick_wand[i] = 0;
    }
  }

  MagickWandTerminus();
  InitMagickWandLib = 0;

exit_exit:

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (rval);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

//
// join wand with default wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".join"
Ent *
ent_im_wand_join (int nargs, Datum args[])
{
  Ent *e1=0, *rent;
  int rval=0, j;

  if (nargs != 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Join or add images from another Magick Wand to default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(another_wand),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'another_wand' contains images that are inserted into default Magick Wand depending on the iterator.\n");
    rerror ( THIS_SOLVER ": requires 1 argument");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": First argument 'wand' must be an integer!");
  j = (int) class_double(e1) - 1;

  if (j<0 || j>=MAX_NUMBER_MAGICK_WANDS || default_wand_idx==j )
    fprintf(stderr,THIS_SOLVER  ": Index out of bounds!\n");

  if (magick_wand[j])
  {
    if(MagickAddImage(magick_wand[default_wand_idx],magick_wand[j]) == MagickFalse)
    {
      WriteWandException(magick_wand[default_wand_idx]);
      rval = 1;
    }
  }

  ent_Clean(e1);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (rval);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

//
// join wand with default wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".append"
Ent *
ent_im_wand_append (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  int j;
  char *c=0;
  MagickBooleanType stack=MagickFalse;

  if (nargs != 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Append images in default Magick Wand from active image onwards to a new Magick Wand,\n");
    fprintf (stdout,
             THIS_SOLVER ": and makes the new Magick Wand default.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(new_wand,stack),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'new_wand' contains images that are stacked depending on\n");
    fprintf (stdout,
             THIS_SOLVER ":   'stack' left-to-right (stack='l'), or 'top-to-bottom' (stack='t').\n");
    fprintf (stdout,
             THIS_SOLVER ": The stacking starts from active image in default wand.\n");
    rerror ( THIS_SOLVER ": requires 2 arguments");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": First argument 'new_wand' must be integer scalar!");
  j = (int) class_double(e1) - 1;

  if (j<0 || j>=MAX_NUMBER_MAGICK_WANDS || default_wand_idx==j )
    fprintf(stderr,THIS_SOLVER  ": Index out of bounds!\n");

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) != MATRIX_DENSE_STRING)
    rerror(THIS_SOLVER ": Second argument 'stack' must be integer scalar!");
  c = class_char_pointer(e2);

  if (c)
  {
    if (*c == 't' || *c == 'T')
      stack = MagickTrue;
  }

  // if wand j exists destroy it
  if (magick_wand[j])
  {
    magick_wand[j] = DestroyMagickWand(magick_wand[j]);
    magick_wand[j] = 0;
  }

  magick_wand[j] = MagickAppendImages(magick_wand[default_wand_idx], stack);
  default_wand_idx = j;

  ent_Clean(e1);
  ent_Clean(e2);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

//
// combine images in default wand into a single image of new (default) wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".combine"
Ent *
ent_im_wand_combine (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  int j;
  char *c=0;
  int ar=0, ac=0, ag=0, am=0, ab=0, ay=0, ao=0, aB=0, aM=0;

  if (nargs != 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Combine images from default Magick Wand to a new Magick Wand,\n");
    fprintf (stdout,
             THIS_SOLVER ": and make the new Magick Wand default.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(new_wand],channel]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'new_wand' contains images that are stacked depending on\n");
    fprintf (stdout,
             THIS_SOLVER ":   'channel' string containing first names of the channels.\n");
    fprintf (stdout,
             THIS_SOLVER ": The stacking starts from active image in default wand.\n");
    rerror ( THIS_SOLVER ": requires 1 or 2 arguments");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": First argument 'new_wand' must be integer scalar!");
  j = (int) class_double(e1) - 1;

  if (j<0 || j>=MAX_NUMBER_MAGICK_WANDS || default_wand_idx==j )
    fprintf(stderr,THIS_SOLVER  ": Index out of bounds!\n");

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": Second argument 'channel' must be integer scalar!");
    c = class_char_pointer(e2);
  }

  if (c)
  {
    //
    // process channels if necessary
    //
    ar = (strstr(c, "r") != 0);  // RedChannel
    ac = (strstr(c, "c") != 0);  // CyanChannel
    ag = (strstr(c, "g") != 0);  // GreenChannel
    am = (strstr(c, "m") != 0);  // MagentaChannel
    ab = (strstr(c, "b") != 0);  // BlueChannel
    ay = (strstr(c, "y") != 0);  // YellowChannel
    ao = (strstr(c, "o") != 0 || strstr(c, "a") != 0);  // OpacityChannel
    aB = (strstr(c, "B") != 0 || strstr(c, "k") != 0);  // BlackChannel
    aM = (strstr(c, "M") != 0 || strstr(c, "t") != 0);  // MatteChannel
  }

  // if wand j exists destroy it
  if (magick_wand[j])
  {
    magick_wand[j] = DestroyMagickWand(magick_wand[j]);
    magick_wand[j] = 0;
  }

  if (ar)
    magick_wand[j] = MagickCombineImages(magick_wand[default_wand_idx], RedChannel);
  else if (ac)
    magick_wand[j] = MagickCombineImages(magick_wand[default_wand_idx], CyanChannel);
  else if (ag)
    magick_wand[j] = MagickCombineImages(magick_wand[default_wand_idx], GreenChannel);
  else if (am)
    magick_wand[j] = MagickCombineImages(magick_wand[default_wand_idx], MagentaChannel);
  else if (ab)
    magick_wand[j] = MagickCombineImages(magick_wand[default_wand_idx], BlueChannel);
  else if (ay)
    magick_wand[j] = MagickCombineImages(magick_wand[default_wand_idx], YellowChannel);
  else if (ao)
    magick_wand[j] = MagickCombineImages(magick_wand[default_wand_idx], OpacityChannel);
  else if (aB)
    magick_wand[j] = MagickCombineImages(magick_wand[default_wand_idx], BlackChannel);
  else if (aM)
    magick_wand[j] = MagickCombineImages(magick_wand[default_wand_idx], MatteChannel);
  else if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
    magick_wand[j] = MagickCombineImages(magick_wand[default_wand_idx], DefaultChannels);

  default_wand_idx = j;

  ent_Clean(e1);
  ent_Clean(e2);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

//
// combine images in default wand into a single image of new (default) wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".diff"
Ent *
ent_im_wand_diff (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  int j,n;
  char *c=0;
  int ar=0, ac=0, ag=0, am=0, ab=0, ay=0, ao=0, aB=0, aM=0;
  MetricType metric=AbsoluteErrorMetric;;

  double distortion;

  if (nargs != 3 && nargs != 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Compare images from default Magick Wand to another Magick Wand,\n");
    fprintf (stdout,
             THIS_SOLVER ": and put result in new Magick Wand that is also default.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(another_wand,metric,[channel]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'another_wand' contains images that are stacked depending on\n");
    fprintf (stdout,
             THIS_SOLVER ":   'channel' string containing first names of the channels.\n");
    fprintf (stdout,
             THIS_SOLVER ": The stacking starts from active image in default wand.\n");
    rerror ( THIS_SOLVER ": requires 1 or 2 arguments");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": First argument 'another_wand' must be integer scalar!");
  j = (int) class_double(e1) - 1;

  if (j<0 || j>=MAX_NUMBER_MAGICK_WANDS || default_wand_idx==j || !magick_wand[j])
  {
    fprintf(stderr,THIS_SOLVER  ": Index out of bounds!\n");
    rerror(THIS_SOLVER ": First argument 'another_wand' must be integer scalar!");
  }

  for (n=0; n<MAX_NUMBER_MAGICK_WANDS; n++)
  {
    if (magick_wand[n])
      continue;
    else
      break;
  }
  if (n==MAX_NUMBER_MAGICK_WANDS)
  {
    fprintf(stderr,THIS_SOLVER  ": All the wands are busy serving other customers! Get rid of some!\n");
    rerror(THIS_SOLVER ": Cannot continue!");
  }

  if (nargs>=2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": Second argument 'metric' must be string!");
    c = class_char_pointer(e2);
    if (!c)
      rerror(THIS_SOLVER ": Second argument 'metric' must be string!");

    if (strstr(c,"abs"))
      metric = AbsoluteErrorMetric;
    else if (strstr(c,"meanabs"))
      metric = MeanAbsoluteErrorMetric;
    else if (strstr(c,"meanpp"))
      metric = MeanErrorPerPixelMetric;
    else if (strstr(c,"square"))
      metric = MeanSquaredErrorMetric;
    else if (strstr(c,"peakabs"))
      metric = PeakAbsoluteErrorMetric;
    else if (strstr(c,"peaksnr"))
      metric = PeakSignalToNoiseRatioMetric;
    else if (strstr(c,"rms"))
      metric = RootMeanSquaredErrorMetric;
    else if (strstr(c,"normcc"))
      metric = NormalizedCrossCorrelationErrorMetric;
    else if (strstr(c,"fuzz"))
      metric = FuzzErrorMetric;
    else
      rerror(THIS_SOLVER ": Second argument 'metric' out of range!");

    c = 0;
  }

  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3) != MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": Second argument 'channel' must be integer scalar!");
    c = class_char_pointer(e3);
  }
  if (c)
  {
    //
    // process channels if necessary
    //
    ar = (strstr(c, "r") != 0);  // RedChannel
    ac = (strstr(c, "c") != 0);  // CyanChannel
    ag = (strstr(c, "g") != 0);  // GreenChannel
    am = (strstr(c, "m") != 0);  // MagentaChannel
    ab = (strstr(c, "b") != 0);  // BlueChannel
    ay = (strstr(c, "y") != 0);  // YellowChannel
    ao = (strstr(c, "o") != 0) || (strstr(c, "a") != 0);  // OpacityChannel
    aB = (strstr(c, "B") != 0) || (strstr(c, "k") != 0);  // BlackChannel
    aM = (strstr(c, "M") != 0) || (strstr(c, "t") != 0);  // MatteChannel
  }

  // Do the work now
  if (ar)
    magick_wand[n] = MagickCompareImageChannels(magick_wand[default_wand_idx], magick_wand[j],
                                                RedChannel, metric, &distortion);
  if (ac)
    magick_wand[n] = MagickCompareImageChannels(magick_wand[default_wand_idx], magick_wand[j],
                                                CyanChannel, metric, &distortion);
  if (ag)
    magick_wand[n] = MagickCompareImageChannels(magick_wand[default_wand_idx], magick_wand[j],
                                                GreenChannel, metric, &distortion);
  if (am)
    magick_wand[n] = MagickCompareImageChannels(magick_wand[default_wand_idx], magick_wand[j],
                                                MagentaChannel, metric, &distortion);
  if (ab)
    magick_wand[n] = MagickCompareImageChannels(magick_wand[default_wand_idx], magick_wand[j],
                                                BlueChannel, metric, &distortion);
  if (ay)
    magick_wand[n] = MagickCompareImageChannels(magick_wand[default_wand_idx], magick_wand[j],
                                                YellowChannel, metric, &distortion);
  if (ao)
    magick_wand[n] = MagickCompareImageChannels(magick_wand[default_wand_idx], magick_wand[j],
                                                OpacityChannel, metric, &distortion);
  if (aB)
    magick_wand[n] = MagickCompareImageChannels(magick_wand[default_wand_idx], magick_wand[j],
                                                BlackChannel, metric, &distortion);
  if (aM)
    magick_wand[n] = MagickCompareImageChannels(magick_wand[default_wand_idx], magick_wand[j],
                                                MatteChannel, metric, &distortion);
  if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
  {
    magick_wand[n] = MagickCompareImages(magick_wand[default_wand_idx], magick_wand[j],
                                         metric, &distortion);
  }

  if (magick_wand[n])
  {
    default_wand_idx = n;

    // add distortion to the artifacts of default image
    sprintf(default_text, "%g", distortion);
    MagickSetImageArtifact(magick_wand[default_wand_idx], "distortion", default_text);
  }
  else
  {
    fprintf(stderr,THIS_SOLVER ": Failed to create new wand!\n");
    rerror (THIS_SOLVER ": Terrible Internal Error. Cannot continue!\n");
  }

  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}


//
// combine images in default wand into a single image of new (default) wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".composite"
Ent *
ent_im_wand_composite (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *x=0;
  int n;
  char *c=0;
  int ar=0, ac=0, ag=0, am=0, ab=0, ay=0, ao=0, aB=0, aM=0;

  CompositeOperator compose=ModulusAddCompositeOp;

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs > 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Composite images from default Magick Wand to a new Magick Wand,\n");
    fprintf (stdout,
             THIS_SOLVER ": and make the new Magick Wand default.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(op,offset,[channel]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'op' is string containing the composition operator,\n");
    fprintf (stdout,
             THIS_SOLVER ":   'offset' [x,y] is two-component vector, while\n");
    fprintf (stdout,
             THIS_SOLVER ":   'channel' string containing first names of the channels.\n");
    fprintf (stdout,
             THIS_SOLVER ": The stacking starts from active image in default wand.\n");
    rerror ( THIS_SOLVER ": requires no or 1 argument");
  }

  if (nargs >= 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) != MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": First argument 'op' must be string!");
    c = class_char_pointer(e1);

    if (!strcmp(c,"mod"))
      compose = ModulusAddCompositeOp;
    else if (!strcmp(c,"atop"))
      compose = AtopCompositeOp;
    else if (!strcmp(c,"blend"))
      compose = BlendCompositeOp;
    else if (!strcmp(c,"bumpmap"))
      compose = BumpmapCompositeOp;
    else if (!strcmp(c,"changemask"))
      compose = ChangeMaskCompositeOp;
    else if (!strcmp(c,"clear"))
      compose = ClearCompositeOp;
    else if (!strcmp(c,"burn"))
      compose = ColorBurnCompositeOp;
    else if (!strcmp(c,"dodge"))
      compose = ColorDodgeCompositeOp;
    else if (!strcmp(c,"colorize"))
      compose = ColorizeCompositeOp;
    else if (!strcmp(c,"copyblack"))
      compose = CopyBlackCompositeOp;
    else if (!strcmp(c,"copyblue"))
      compose = CopyBlueCompositeOp;
    else if (!strcmp(c,"copy"))
      compose = CopyCompositeOp;
    else if (!strcmp(c,"copycyan"))
      compose = CopyCyanCompositeOp;
    else if (!strcmp(c,"copygreen"))
      compose = CopyGreenCompositeOp;
    else if (!strcmp(c,"copymag"))
      compose = CopyMagentaCompositeOp;
    else if (!strcmp(c,"copyopac"))
      compose = CopyOpacityCompositeOp;
    else if (!strcmp(c,"copyred"))
      compose = CopyRedCompositeOp;
    else if (!strcmp(c,"copyyell"))
      compose = CopyYellowCompositeOp;
    else if (!strcmp(c,"darken"))
      compose = DarkenCompositeOp;
    else if (!strcmp(c,"dstatop"))
      compose = DstAtopCompositeOp;
    else if (!strcmp(c,"dst"))
      compose = DstCompositeOp;
    else if (!strcmp(c,"dstin"))
      compose = DstInCompositeOp;
    else if (!strcmp(c,"dstout"))
      compose = DstOutCompositeOp;
    else if (!strcmp(c,"dstover"))
      compose = DstOverCompositeOp;
    else if (!strcmp(c,"diff"))
      compose = DifferenceCompositeOp;
    else if (!strcmp(c,"displace"))
      compose = DisplaceCompositeOp;
    else if (!strcmp(c,"dissolv"))
      compose = DissolveCompositeOp;
    else if (!strcmp(c,"excl"))
      compose = ExclusionCompositeOp;
    else if (!strcmp(c,"hardlight"))
      compose = HardLightCompositeOp;
    else if (!strcmp(c,"huecomp"))
      compose = HueCompositeOp;
    else if (!strcmp(c,"in"))
      compose = InCompositeOp;
    else if (!strcmp(c,"lighten"))
      compose = LightenCompositeOp;
    else if (!strcmp(c,"linlight"))
      compose = LinearLightCompositeOp;
    else if (!strcmp(c,"luminize"))
      compose = LuminizeCompositeOp;
    else if (!strcmp(c,"minusdst"))
      compose = MinusDstCompositeOp;
    else if (!strcmp(c,"modulate"))
      compose = ModulateCompositeOp;
    else if (!strcmp(c,"multiply"))
      compose = MultiplyCompositeOp;
    else if (!strcmp(c,"out"))
      compose = OutCompositeOp;
    else if (!strcmp(c,"over"))
      compose = OverCompositeOp;
    else if (!strcmp(c,"overlay"))
      compose = OverlayCompositeOp;
    else if (!strcmp(c,"plus"))
      compose = PlusCompositeOp;
    else if (!strcmp(c,"replace"))
      compose = ReplaceCompositeOp;
    else if (!strcmp(c,"saturate"))
      compose = SaturateCompositeOp;
    else if (!strcmp(c,"screen"))
      compose = ScreenCompositeOp;
    else if (!strcmp(c,"softlight"))
      compose = SoftLightCompositeOp;
    else if (!strcmp(c,"srcatop"))
      compose = SrcAtopCompositeOp;
    else if (!strcmp(c,"src"))
      compose = SrcCompositeOp;
    else if (!strcmp(c,"srcin"))
      compose = SrcInCompositeOp;
    else if (!strcmp(c,"srcout"))
      compose = SrcOutCompositeOp;
    else if (!strcmp(c,"srcover"))
      compose = SrcOverCompositeOp;
    else if (!strcmp(c,"modsub"))
      compose = ModulusSubtractCompositeOp;
    else if (!strcmp(c,"threshold"))
      compose = ThresholdCompositeOp;
    else if (!strcmp(c,"xor"))
      compose = XorCompositeOp;
    else if (!strcmp(c,"divide"))
      compose = DivideDstCompositeOp;
    else if (!strcmp(c,"distort"))
      compose = DistortCompositeOp;
    else if (!strcmp(c,"blur"))
      compose = BlurCompositeOp;
    else if (!strcmp(c,"pegtop"))
      compose = PegtopLightCompositeOp;
    else if (!strcmp(c,"vividlight"))
      compose = VividLightCompositeOp;
    else if (!strcmp(c,"pinlight"))
      compose = PinLightCompositeOp;
    else if (!strcmp(c,"lindodge"))
      compose = LinearDodgeCompositeOp;
    else if (!strcmp(c,"linburn"))
      compose = LinearBurnCompositeOp;
    else if (!strcmp(c,"math"))
      compose = MathematicsCompositeOp;
    else if (!strcmp(c,"divsrc"))
      compose = DivideSrcCompositeOp;
    else if (!strcmp(c,"minsrc"))
      compose = MinusSrcCompositeOp;
    else if (!strcmp(c,"darkenint"))
      compose = DarkenIntensityCompositeOp;
    else if (!strcmp(c,"lightenint"))
      compose = LightenIntensityCompositeOp;
    else
      rerror(THIS_SOLVER ": I cannot recognize first argument 'op'!");

    c = 0;
  }

  if (nargs >= 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": Second argument 'offset' must be two-column vector!");

    x = ent_data(e2);
  }

  if (nargs >= 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3) != MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": Third argument 'channel' must be integer scalar!");
    c = class_char_pointer(e3);
  }
  if (c)
  {
    //
    // process channels if necessary
    //
    ar = (strstr(c, "r") != 0);  // RedChannel
    ac = (strstr(c, "c") != 0);  // CyanChannel
    ag = (strstr(c, "g") != 0);  // GreenChannel
    am = (strstr(c, "m") != 0);  // MagentaChannel
    ab = (strstr(c, "b") != 0);  // BlueChannel
    ay = (strstr(c, "y") != 0);  // YellowChannel
    ao = (strstr(c, "o") != 0);  // OpacityChannel
    aB = (strstr(c, "B") != 0);  // BlackChannel
    aM = (strstr(c, "M") != 0);  // MatteChannel
  }

  for (n=0; n<MAX_NUMBER_MAGICK_WANDS; n++)
  {
    if (magick_wand[n])
      continue;
    else
      break;
  }
  if (n==MAX_NUMBER_MAGICK_WANDS)
  {
    fprintf(stderr,THIS_SOLVER  ": All the wands are busy serving other customers! Get rid of some!\n");
    rerror(THIS_SOLVER ": Cannot continue!");
  }

  magick_wand[n] = NewMagickWand();

  if (ar)
    MagickCompositeImageChannel(magick_wand[n], RedChannel, magick_wand[default_wand_idx], compose,
                                (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));
  else if (ac)
    MagickCompositeImageChannel(magick_wand[n], CyanChannel, magick_wand[default_wand_idx], compose,
                                (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));
  else if (ag)
    MagickCompositeImageChannel(magick_wand[n], GreenChannel, magick_wand[default_wand_idx], compose,
                                (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));
  else if (am)
    MagickCompositeImageChannel(magick_wand[n], MagentaChannel, magick_wand[default_wand_idx], compose,
                                (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));
  else if (ab)
    MagickCompositeImageChannel(magick_wand[n], BlueChannel, magick_wand[default_wand_idx], compose,
                                (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));
  else if (ay)
    MagickCompositeImageChannel(magick_wand[n], YellowChannel, magick_wand[default_wand_idx], compose,
                                (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));
  else if (ao)
    MagickCompositeImageChannel(magick_wand[n], OpacityChannel, magick_wand[default_wand_idx], compose,
                                (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));
  else if (aB)
    MagickCompositeImageChannel(magick_wand[n], BlackChannel, magick_wand[default_wand_idx], compose,
                                (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));
  else if (aM)
    MagickCompositeImageChannel(magick_wand[n], MatteChannel, magick_wand[default_wand_idx], compose,
                                (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));
  else if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
    MagickCompositeImage(magick_wand[n], magick_wand[default_wand_idx], compose,
                         (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));

  default_wand_idx = n;

  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}


//
// deconstruct default wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".deconstruct"
Ent *
ent_im_wand_deconstruct (int nargs, Datum args[])
{
  Ent *rent;

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs != 0)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Replaces images in default Magick Wand with their deconstructed version,\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(),\n");
    rerror ( THIS_SOLVER ": requires no arguments");
  }

  MagickDeconstructImages(magick_wand[default_wand_idx]);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

//
// combine images in default wand into a single image of new (default) wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".distort"
Ent *
ent_im_wand_distort (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *x=0;
  char *c=0;
  size_t n=0;

  DistortImageMethod method=AffineDistortion;
  MagickBooleanType bestfit = MagickFalse;

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs > 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Distort active image from default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(distortion,parameters,[bestfit]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'distortion' is string containing the operation,\n");
    fprintf (stdout,
             THIS_SOLVER ":   'parameters' is real vector, while\n");
    fprintf (stdout,
             THIS_SOLVER ":   'bestfit' to be used if resizing.\n");
    rerror ( THIS_SOLVER ": requires 2 or 3 arguments");
  }

  if (nargs >= 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) != MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": First argument 'distortion' must be string!");
    c = class_char_pointer(e1);

    if (!strcmp(c,"affine"))
      method = AffineDistortion;
    else if (!strcmp(c,"affineproj"))
      method = AffineProjectionDistortion;
    else if (!strcmp(c,"scalerottran"))
      method = ScaleRotateTranslateDistortion;
    else if (!strcmp(c,"perspectivedist"))
      method = PerspectiveDistortion;
    else if (!strcmp(c,"perspectiveproj"))
      method = PerspectiveProjectionDistortion;
    else if (!strcmp(c,"bilinearforward"))
      method = BilinearForwardDistortion;
    else if (!strcmp(c,"bilinear"))
      method = BilinearDistortion;
    else if (!strcmp(c,"bilinearrev"))
      method = BilinearReverseDistortion;
    else if (!strcmp(c,"poly"))
      method = PolynomialDistortion;
    else if (!strcmp(c,"arc"))
      method = ArcDistortion;
    else if (!strcmp(c,"polar"))
      method = PolarDistortion;
    else if (!strcmp(c,"depolar"))
      method = DePolarDistortion;
    else if (!strcmp(c,"cyl"))
      method = Cylinder2PlaneDistortion;
    else if (!strcmp(c,"plane"))
      method = Plane2CylinderDistortion;
    else if (!strcmp(c,"barrel"))
      method = BarrelDistortion;
    else if (!strcmp(c,"barrelinv"))
      method = BarrelInverseDistortion;
    else if (!strcmp(c,"shep"))
      method = ShepardsDistortion;
    else if (!strcmp(c,"resize"))
      method = ResizeDistortion;
    else if (!strcmp(c,"sent"))
      method = SentinelDistortion ;
    else
      rerror(THIS_SOLVER ": I cannot recognize first argument 'distortion'!");

    c = 0;
  }

  if (nargs >= 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": Second argument 'parameters' must be real vector!");
    x = ent_data(e2);

    n = x->nrow * x->ncol;
  }

  if (nargs >= 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3) != MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": Third argument 'bestfit' must be 0 or 1!");
    if (class_double(e3)>0)
      bestfit = MagickTrue;
  }

  MagickDistortImage(magick_wand[default_wand_idx], method, n, MDRPTR(x), bestfit);

  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}


//
// evaluate operator on image
//
// #undef THIS_SOLVER
// #define THIS_SOLVER RLAB_NAME_EMBED_IM ".eval"
// Ent *
// ent_im_wand_eval (int nargs, Datum args[])
// {
//   Ent *e1=0, *e2=0, *e3=0, *rent;
//   double value=0;
//   char *c=0;
//
//   int ar=0, ac=0, ag=0, am=0, ab=0, ay=0, ao=0, aB=0, aM=0;
//
//   MagickEvaluateOperator op=UndefinedEvaluateOperator;
//
//   if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
//     rerror(THIS_SOLVER  ": Help! Library not initialized!\n");
//
//   if (nargs == 0)
//   {
//     fprintf (stdout,
//              THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
//     fprintf (stdout,
//              THIS_SOLVER ": Evaluate operator on active image from default Magick Wand.\n");
//     fprintf (stdout,
//              THIS_SOLVER ": Format:\n");
//     fprintf (stdout,
//              THIS_SOLVER ":   " THIS_SOLVER "(operator,value,channel),\n");
//     fprintf (stdout,
//              THIS_SOLVER ": where\n");
//     fprintf (stdout,
//              THIS_SOLVER ":   'operator' is string containing the operation,\n");
//     fprintf (stdout,
//              THIS_SOLVER ":   'value' is real scalar parameter of the operator, while\n");
//     fprintf (stdout,
//              THIS_SOLVER ":   'channel' is the color channel to which the operator is applied.\n");
//     rerror ( THIS_SOLVER ": requires 1,2 or 3 arguments");
//   }
//
//   if (nargs >= 1)
//   {
//     e1 = bltin_get_ent (args[0]);
//     if (ent_type(e1) != MATRIX_DENSE_STRING)
//       rerror(THIS_SOLVER ": First argument 'operator' must be string!");
//     c = class_char_pointer(e1);
//
//     if (!strcmp(c,"+"))
//       op = AddEvaluateOperator;
//     else if (!strcmp(c,"and"))
//       op = AndEvaluateOperator;
//     else if (!strcmp(c,"/") || !strcmp(c,"divide"))
//       op = DivideEvaluateOperator;
//     else if (!strcmp(c,"<<") || !strcmp(c,"shiftleft"))
//       op = LeftShiftEvaluateOperator;
//     else if (!strcmp(c,"max"))
//       op = MaxEvaluateOperator;
//     else if (!strcmp(c,"min"))
//       op = MinEvaluateOperator;
//     else if (!strcmp(c,"*") || !strcmp(c,"multiply"))
//       op = MultiplyEvaluateOperator;
//     else if (!strcmp(c,"or"))
//       op = OrEvaluateOperator;
//     else if (!strcmp(c,">>") || !strcmp(c,"shiftright"))
//       op = RightShiftEvaluateOperator;
//     else if (!strcmp(c,"set"))
//       op = SetEvaluateOperator;
//     else if (!strcmp(c,"-") || !strcmp(c,"subtract") || !strcmp(c,"minus"))
//       op = SubtractEvaluateOperator;
//     else if (!strcmp(c,"xor"))
//       op = XorEvaluateOperator;
//     else if (!strcmp(c,"pow"))
//       op = PowEvaluateOperator;
//     else if (!strcmp(c,"log"))
//       op = LogEvaluateOperator;
//     else if (!strcmp(c,"threshold"))
//       op = ThresholdEvaluateOperator;
//     else if (!strcmp(c,"thresholdblack"))
//       op = ThresholdBlackEvaluateOperator;
//     else if (!strcmp(c,"thresholdwhite"))
//       op = ThresholdWhiteEvaluateOperator;
//     else if (!strcmp(c,"gaussian") || !strcmp(c,"noisegaussian"))
//       op = GaussianNoiseEvaluateOperator;
//     else if (!strcmp(c,"impulse") || !strcmp(c,"noiseimpulse"))
//       op = ImpulseNoiseEvaluateOperator;
//     else if (!strcmp(c,"laplace") || !strcmp(c,"noiselaplace"))
//       op = LaplacianNoiseEvaluateOperator;
//     else if (!strcmp(c,"multiplicative") || !strcmp(c,"noisemultiplicative"))
//       op = MultiplicativeNoiseEvaluateOperator;
//     else if (!strcmp(c,"poisson") || !strcmp(c,"noisepoisson"))
//       op = PoissonNoiseEvaluateOperator;
//     else if (!strcmp(c,"uniform") || !strcmp(c,"noiseuniform"))
//       op = UniformNoiseEvaluateOperator;
//     else if (!strcmp(c,"cos"))
//       op = CosineEvaluateOperator;
//     else if (!strcmp(c,"sin"))
//       op = SineEvaluateOperator;
//     else if (!strcmp(c,"+mod") || !strcmp(c,"addmodulus"))
//       op = AddModulusEvaluateOperator;
//     else if (!strcmp(c,"mean"))
//       op = MeanEvaluateOperator;
//     else if (!strcmp(c,"abs"))
//       op = AbsEvaluateOperator;
//     else if (!strcmp(c,"exp"))
//       op = ExponentialEvaluateOperator;
//     else if (!strcmp(c,"median"))
//       op = MedianEvaluateOperator;
// //     else if (!strcmp(c,"sum"))
// //       op = SumEvaluateOperator;
//     else
//       rerror(THIS_SOLVER ": I cannot recognize first argument 'operator'!");
//
//     c = 0;
//   }
//
//   if (nargs >= 2)
//   {
//     e2 = bltin_get_ent (args[1]);
//     if (ent_type(e2) != MATRIX_DENSE_REAL)
//       rerror(THIS_SOLVER ": Second argument 'value' must be real scalar!");
//     value = class_double(e2);
//   }
//
//   if (nargs >= 3)
//   {
//     e3 = bltin_get_ent (args[2]);
//     if (ent_type(e3) != MATRIX_DENSE_STRING)
//       rerror(THIS_SOLVER ": Third argument 'channel' must be integer scalar!");
//     c = class_char_pointer(e3);
//   }
//   if (c)
//   {
//     //
//     // process channels if necessary
//     //
//     ar = (strstr(c, "r") != 0);  // RedChannel
//     ac = (strstr(c, "c") != 0);  // CyanChannel
//     ag = (strstr(c, "g") != 0);  // GreenChannel
//     am = (strstr(c, "m") != 0);  // MagentaChannel
//     ab = (strstr(c, "b") != 0);  // BlueChannel
//     ay = (strstr(c, "y") != 0);  // YellowChannel
//     ao = (strstr(c, "o") != 0);  // OpacityChannel
//     aB = (strstr(c, "B") != 0);  // BlackChannel
//     aM = (strstr(c, "M") != 0);  // MatteChannel
//   }
//
//   if (ar)
//     MagickEvaluateImageChannel(magick_wand[default_wand_idx], RedChannel, op, value);
//   else if (ac)
//     MagickEvaluateImageChannel(magick_wand[default_wand_idx], CyanChannel, op, value);
//   else if (ag)
//     MagickEvaluateImageChannel(magick_wand[default_wand_idx], GreenChannel, op, value);
//   else if (am)
//     MagickEvaluateImageChannel(magick_wand[default_wand_idx], MagentaChannel, op, value);
//   else if (ab)
//     MagickEvaluateImageChannel(magick_wand[default_wand_idx], BlueChannel, op, value);
//   else if (ay)
//     MagickEvaluateImageChannel(magick_wand[default_wand_idx], YellowChannel, op, value);
//   else if (ao)
//     MagickEvaluateImageChannel(magick_wand[default_wand_idx], OpacityChannel, op, value);
//   else if (aB)
//     MagickEvaluateImageChannel(magick_wand[default_wand_idx], BlackChannel, op, value);
//   else if (aM)
//     MagickEvaluateImageChannel(magick_wand[default_wand_idx], MatteChannel, op, value);
//   else if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
//     MagickEvaluateImage(magick_wand[default_wand_idx], op, value);
//
//   ent_Clean(e1);
//   ent_Clean(e2);
//   ent_Clean(e3);
//
//   rent = ent_Create ();
//   ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
//   ent_type (rent) = MATRIX_DENSE_REAL;
//   return (rent);
// }

//
// combine image with an image containing color lookup table
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".clut"
Ent *
ent_im_wand_clut (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  int j;
  char *c=0;
  int ar=0, ac=0, ag=0, am=0, ab=0, ay=0, ao=0, aB=0, aM=0;

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs != 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Convert active image in default Magick Wand using color lookup table from another Magick Wand,\n");
    fprintf (stdout,
             THIS_SOLVER ": and make the new Magick Wand default.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(clut_wand[,channel]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'clut_wand' contains color lookup table, while\n");
    fprintf (stdout,
             THIS_SOLVER ":   'channel' string containing first names of the channels.\n");
    fprintf (stdout,
             THIS_SOLVER ": The stacking starts from active image in default wand.\n");
    rerror ( THIS_SOLVER ": requires 1 or 2 arguments");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": First argument 'clut_wand' must be integer scalar!");
  j = (int) class_double(e1) - 1;

  if (j<0 || j>=MAX_NUMBER_MAGICK_WANDS || default_wand_idx==j || !magick_wand[j])
  {
    fprintf(stderr,THIS_SOLVER  ": Index out of bounds!\n");
    rerror(THIS_SOLVER ": First argument 'clut_wand' must be integer scalar!");
  }

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": Second argument 'channel' must be integer scalar!");
    c = class_char_pointer(e2);
  }

  if (c)
  {
    //
    // process channels if necessary
    //
    ar = (strstr(c, "r") != 0);  // RedChannel
    ag = (strstr(c, "g") != 0);  // GreenChannel
    ab = (strstr(c, "b") != 0);  // BlueChannel
    ao = ((strstr(c, "o") != 0) || (strstr(c, "a") != 0));  // OpacityChannel or alpha
    ac = (strstr(c, "c") != 0);  // CyanChannel
    am = (strstr(c, "m") != 0);  // MagentaChannel
    ay = (strstr(c, "y") != 0);  // YellowChannel
    aB = ((strstr(c, "k") != 0)||(strstr(c, "B") != 0));  // BlackChannel
    aM = ((strstr(c, "M") != 0)||(strstr(c, "t") != 0));  // MatteChannel
  }

  if (ar)
    MagickClutImageChannel(magick_wand[default_wand_idx], RedChannel, magick_wand[j]);
  else if (ac)
    MagickClutImageChannel(magick_wand[default_wand_idx], CyanChannel, magick_wand[j]);
  else if (ag)
    MagickClutImageChannel(magick_wand[default_wand_idx], GreenChannel, magick_wand[j]);
  else if (am)
    MagickClutImageChannel(magick_wand[default_wand_idx], MagentaChannel, magick_wand[j]);
  else if (ab)
    MagickClutImageChannel(magick_wand[default_wand_idx], BlueChannel, magick_wand[j]);
  else if (ay)
    MagickClutImageChannel(magick_wand[default_wand_idx], YellowChannel, magick_wand[j]);
  else if (ao)
    MagickClutImageChannel(magick_wand[default_wand_idx], OpacityChannel, magick_wand[j]);
  else if (aB)
    MagickClutImageChannel(magick_wand[default_wand_idx], BlackChannel, magick_wand[j]);
  else if (aM)
    MagickClutImageChannel(magick_wand[default_wand_idx], MatteChannel, magick_wand[j]);
  else if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
    MagickClutImage(magick_wand[default_wand_idx], magick_wand[j]);

  ent_Clean(e1);
  ent_Clean(e2);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

// #undef THIS_SOLVER
// #define THIS_SOLVER RLAB_NAME_EMBED_IM ".constitute"
// Ent *
// ent_im_wand_image_constitute(int nargs, Datum args[])
// {
//   Ent *e1=0, *en=0, *rent;
//   ListNode *node, *node2;
//   MDR *pixel_data[4]={0};
//
//   char channel[2];
//   unsigned char pixel_storage=0, is_rgba=0, is_cmyk=0, is_gray=0, is_pad=0;
//
//   double ddummy;
//   int npf=0, i, j, k;
//
//   // image info:
//   int width=-1,height=-1;
//   int lwidth=-1,lheight=-1;
//   double qrange=RLAB_MAGICK_QUANTUM_SIZE;
//   char *map=0;
//   unsigned char *pixels_8=0;
//   short unsigned int *pixels_16=0;
//
//   if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
//     rerror(THIS_SOLVER  ": Help! Library not initialized!\n");
//
//   if (nargs != 1)
//   {
//     fprintf (stdout,
//              THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
//     fprintf (stdout,
//              THIS_SOLVER ": Constitute new image in default Magick Wand from the pixel data provided,\n");
//     fprintf (stdout,
//              THIS_SOLVER ": Format:\n");
//     fprintf (stdout,
//              THIS_SOLVER ":   " THIS_SOLVER "(pixel_data),\n");
//     fprintf (stdout,
//              THIS_SOLVER ": where\n");
//     fprintf (stdout,
//              THIS_SOLVER ":   'pixel_data' is a list <<pixel;pixel_map[;height][;width][;range]>>,\n");
//     fprintf (stdout,
//              THIS_SOLVER ":   containing mapping and data for all pixels in the image.\n");
//     fprintf (stdout,
//              THIS_SOLVER ": See manual for more info.\n");
//     rerror ( THIS_SOLVER ": requires 1 or 2 arguments");
//   }
//
//   e1 = bltin_get_ent (args[0]);
//   if (ent_type (e1) != BTREE)
//     rerror(THIS_SOLVER ": First argument 'image' must be list <<pixel;pixel_map>>");
//
//   // pixel_map
//   node = btree_FindNode (ent_data (e1), "qrange");
//   if (!node)
//   {
//     ddummy = class_double (var_ent (node));
//     if (ddummy==8 || ddummy==255 || ddummy==1)
//       qrange = 255;
//     else if (ddummy==16 || ddummy==65535 || ddummy==2)
//       qrange = 65535;
//   }
//   // height?
//   node = btree_FindNode (ent_data (e1), "height");
//   if (node != 0)
//   {
//     ddummy = class_double (var_ent (node));
//     if (ddummy > 0 && ddummy==ceil(ddummy))
//       height = (size_t) ddummy;
//   }
//   // width?
//   node = btree_FindNode (ent_data (e1), "width");
//   if (node != 0)
//   {
//     ddummy = class_double (var_ent (node));
//     if (ddummy > 0 && ddummy==ceil(ddummy))
//       width = (size_t) ddummy;
//   }
//   // pixel_map
//   node = btree_FindNode (ent_data (e1), "pixel_map");
//   if (node)
//     map = class_char_pointer (var_ent (node));
//   if (!map)
//   {
//     fprintf(stderr, THIS_SOLVER ": First argument 'image' does not have entry 'pixel_map'\n");
//     rerror(THIS_SOLVER ": Cannot continue!\n");
//   }
//   // pixel is a list containing entries r,g,b,a or c,m,y,k or i or p
//   node = btree_FindNode (ent_data (e1), "pixel");
//   if (!node)
//   {
//     fprintf(stderr, THIS_SOLVER ": First argument 'image' does not have list 'pixel'\n");
//     rerror(THIS_SOLVER ": Cannot continue!\n");
//   }
//
//   en = var_ent (node);
//   if (ent_type (en) != BTREE)
//     rerror(THIS_SOLVER ": Entry 'pixel' must be list <<ch_1;ch_2;..>>");
//
//   npf = strlen(map);
//   for (i=0; i<npf; i++)
//   {
//
//     // map contains initials of channels. we expect these to appear in the list pixel
//     channel[0] = map[i];
//     channel[1] = '\0';
//     node2 = btree_FindNode (ent_data (en), channel);
//     if (!node2)
//     {
//       fprintf(stderr, THIS_SOLVER ": Inconsistency between 'pixel' and 'pixel_map': channel '%s' missing from list 'pixel'!\n",
//               channel);
//       rerror(THIS_SOLVER ": Cannot continue!\n");
//     }
//
//     pixel_data[i] = (MDR *) ent_data(var_ent(node2));
//     if (pixel_data[i]->type==RLAB_TYPE_INT32  && (pixel_storage==0 || pixel_storage=='i'))
//       pixel_storage = 'i';
//     if (pixel_data[i]->type==RLAB_TYPE_DOUBLE && (pixel_storage==0 || pixel_storage=='d'))
//       pixel_storage = 'd';
//
//     if (pixel_storage==0)
//     {
//       fprintf(stderr, THIS_SOLVER ": Inconsistency between 'pixel' and 'pixel_map': channel '%s' missing from list 'pixel'!\n",
//               channel);
//       rerror(THIS_SOLVER ": Cannot continue!\n");
//     }
//     lwidth  = pixel_data[i]->ncol > lwidth  ? pixel_data[i]->ncol : lwidth;
//     lheight = pixel_data[i]->nrow > lheight ? pixel_data[i]->nrow : lheight;
//   }
//
//   // check that mapping is correct
//   if ((!strstr(map,"r") || !strstr(map,"g") || !strstr(map,"b") || !strstr(map,"a")) && !is_cmyk && !is_gray && !is_pad)
//     is_rgba = 1;
//   if ((!strstr(map,"c") || !strstr(map,"m") || !strstr(map,"y") || !strstr(map,"k")) && !is_rgba && !is_gray && !is_pad)
//     is_cmyk = 1;
//   if (!strstr(map,"i") && !is_rgba && !is_cmyk && !is_pad)
//     is_gray = 1;
//   if (!strstr(map,"p") && !is_rgba && !is_cmyk && !is_gray)
//     is_pad = 1;
//   if (is_cmyk + is_rgba + is_gray + is_pad != 1)
//   {
//     fprintf(stderr, THIS_SOLVER ": List 'image' has improper 'pixel_map'\n");
//     rerror(THIS_SOLVER ": Cannot continue!\n");
//   }
//
//   // fix the size
//   if (height == -1)
//     height = lheight;
//   if (width == -1)
//     width = lwidth;
//
//   if (qrange==255)
//   {
//     pixels_8 = (unsigned char *) GC_malloc(height * width * npf * sizeof(unsigned char));
//     for (i=0; i<height; i++)
//     {
//       for (j=0; j<width; j++)
//       {
//         for (k=0; k<npf; k++)
//         {
//           if (pixel_storage == 'd')
//           {
//             pixels_8[i*npf*width + npf*j + k] = (unsigned char)(qrange * Mdr0( pixel_data[k],i,j));
//           }
//           else if (pixel_storage == 'i')
//             pixels_8[i*npf*width + npf*j + k] = (unsigned char) Mdi0(pixel_data[k],i,j) > (unsigned char) qrange ? (unsigned char)qrange : (unsigned char)Mdi0(pixel_data[k],i,j);
//         }
//       }
//     }
//     MagickConstituteImage(magick_wand[default_wand_idx], width, height, map, CharPixel, (void *)pixels_8);
//
//     GC_free(pixels_8);
//   }
//   else if (qrange==65535)
//   {
//     pixels_16 = GC_malloc(height * width * npf * sizeof(short unsigned int));
//
//     for (j=0; j<height; j++)
//     {
//       for (i=0; i<width; i++)
//       {
//         for (k=0; k<npf; k++)
//         {
//           if (pixel_storage == 'd')
//             pixels_16[i*npf*width + npf*j + k] = (unsigned short int)(qrange * Mdr0(pixel_data[k],i,j));
//           else if (pixel_storage == 'i')
//             pixels_16[i*npf*width + npf*j + k] = (unsigned short int) Mdi0(pixel_data[k],i,j) > (unsigned short int) qrange
//                 ? (unsigned short int) qrange : (unsigned short int) Mdi0(pixel_data[k],i,j);
//         }
//       }
//     }
//     MagickConstituteImage(magick_wand[default_wand_idx], (size_t) width, (size_t) height, map, ShortPixel, (void *)pixels_16);
//
//     GC_free(pixels_16);
//   }
//
//   ent_Clean (e1);
//   ent_Clean (en);
//
//   rent = ent_Create ();
//   ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
//   ent_type (rent) = MATRIX_DENSE_REAL;
//   return (rent);
// }


#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".pixels"
Ent * ent_im_wand_image_pixels(int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *en=0, *rent;
  ListNode *node, *node2;
  MDR *pixel_data[5]={0};

  char channel[2];
  unsigned char pixel_storage=0, is_rgba=0, is_cmyka=0, is_gray=0, is_pad=0;

  double ddummy;
  int npf=0, i, j, k;

  // image info:
  int width=-1,height=-1;
  int lwidth=-1,lheight=-1;
  int x=0,y=0;
  double qrange=RLAB_MAGICK_QUANTUM_SIZE;
  char *map=0;
  unsigned char *pixels_8=0;
  short unsigned int *pixels_16=0;

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs != 1 && nargs !=3 && nargs !=4)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Constitute new image in default Magick Wand from the pixel data provided,\n");
    fprintf (stdout,
             THIS_SOLVER ": or extract pixel data from active image in default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ": (1) " THIS_SOLVER "(pixel_data),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'pixel_data' is a list <<pixel;pixel_map[;height][;width][;range]>>,\n");
    fprintf (stdout,
             THIS_SOLVER ":   containing mapping and data for all pixels in the image.\n");
    fprintf (stdout,
             THIS_SOLVER ": (2) pic = " THIS_SOLVER "([y1,y2],[x1,x2],pixel_map[,qrange]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   '[y,w]' is range of pixel columns to be extracted,\n");
    fprintf (stdout,
             THIS_SOLVER ":   '[x,h]' is range of pixel rows to be extracted,\n");
    fprintf (stdout,
             THIS_SOLVER ":   'pixel_map' determines the color map of the extracted pixel matrices, while\n");
    fprintf (stdout,
             THIS_SOLVER ":   'qrange' determines quantum range of the colors.\n");
    fprintf (stdout,
             THIS_SOLVER ": See manual for more info.\n");
    rerror ( THIS_SOLVER ": requires 1 or 2 arguments");
  }

  if (nargs == 1)
  {
    //
    // put pixels into image
    //
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != BTREE)
      rerror(THIS_SOLVER ": First argument 'image' must be list <<pixel;pixel_map>>");

    // pixel_map
    node = btree_FindNode (ent_data (e1), RLAB_NAME_IMG_QRANGE);
    if (node)
    {
      if (ent_type(var_ent(node))==MATRIX_DENSE_REAL)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy==8 || ddummy==255 || ddummy==1)
          qrange = 255;
        else if (ddummy==16 || ddummy==65535 || ddummy==2)
          qrange = 65535;
      }
    }

    // x?
    node = btree_FindNode (ent_data (e1), RLAB_NAME_IMG_X);
    if (node)
    {
      if (ent_type(var_ent(node))==MATRIX_DENSE_REAL)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0 && ddummy==ceil(ddummy))
          x = (int) ddummy - 1;
      }
    }

    // y?
    node = btree_FindNode (ent_data (e1), RLAB_NAME_IMG_Y);
    if (node != 0)
    {
      if (ent_type(var_ent(node))==MATRIX_DENSE_REAL)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0 && ddummy==ceil(ddummy))
          y = (int) ddummy - 1;
      }
    }

    // height?
    node = btree_FindNode (ent_data (e1), RLAB_NAME_IMG_HEIGHT);
    if (node != 0)
    {
      if (ent_type(var_ent(node))==MATRIX_DENSE_REAL)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0 && ddummy==ceil(ddummy))
          height = (int) ddummy;
      }
    }

    // width?
    node = btree_FindNode (ent_data (e1), RLAB_NAME_IMG_WIDTH);
    if (node != 0)
    {
      if (ent_type(var_ent(node))==MATRIX_DENSE_REAL)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0 && ddummy==ceil(ddummy))
          width = (int) ddummy;
      }
    }

    // pixel_map
    node = btree_FindNode (ent_data (e1), RLAB_NAME_IMG_PIXEL_MAP);
    if (node)
    {
      if (ent_type(var_ent(node))==MATRIX_DENSE_STRING)
      {
        map = class_char_pointer (var_ent (node));
      }
    }
    if (!map)
    {
      fprintf(stderr, THIS_SOLVER ": First argument 'image' does not have entry 'pixel_map'\n");
      rerror(THIS_SOLVER ": Cannot continue!\n");
    }
    if (isvalidstring(map)<1 || isvalidstring(map)>5)
    {
      fprintf(stderr, THIS_SOLVER ": Entry 'pixel_map' can be between 1 and 4 character long!\n");
      rerror(THIS_SOLVER ": Cannot continue!\n");
    }


    // pixel is a list containing entries r,g,b,a or c,m,y,k or i or p
    node = btree_FindNode (ent_data (e1), RLAB_NAME_IMG_PIXEL);
    if (!node)
    {
      fprintf(stderr, THIS_SOLVER ": First argument 'image' does not have list 'pixel'\n");
      rerror(THIS_SOLVER ": Cannot continue!\n");
    }

    en = var_ent (node);
    if (ent_type (en) != BTREE)
      rerror(THIS_SOLVER ": Entry 'pixel' must be list <<ch_1;ch_2;..>>");

    npf = strlen(map);
    pixel_storage = 0;
    for (i=0; i<npf; i++)
    {

      // map contains initials of channels. we expect these to appear in the list pixel
      channel[0] = map[i];
      channel[1] = '\0';
      node2 = btree_FindNode (ent_data (en), channel);
      if (!node2)
      {
        fprintf(stderr, THIS_SOLVER ": Inconsistency between 'pixel' and 'pixel_map': channel '%s' missing from list 'pixel'!\n",
                channel);
        rerror(THIS_SOLVER ": Cannot continue!\n");
      }

      pixel_data[i] = (MDR *) ent_data(var_ent(node2));
      if (pixel_data[i]->type==RLAB_TYPE_INT32  && (pixel_storage==0 || pixel_storage=='i'))
        pixel_storage = 'i';
      if (pixel_data[i]->type==RLAB_TYPE_DOUBLE && (pixel_storage==0 || pixel_storage=='d'))
        pixel_storage = 'd';

      if (pixel_storage==0)
      {
        fprintf(stderr, THIS_SOLVER ": Inconsistency between 'pixel' and 'pixel_map': channel '%s' missing from list 'pixel'!\n",
                channel);
        rerror(THIS_SOLVER ": Cannot continue!\n");
      }
      lwidth  = MNC(pixel_data[i]) > lwidth  ? MNC(pixel_data[i]) : lwidth;
      lheight = MNR(pixel_data[i]) > lheight ? MNR(pixel_data[i]) : lheight;
    }


    // check that mapping is correct
    if ((!strstr(map,"r") || !strstr(map,"g") || !strstr(map,"b")) && !is_cmyka && !is_gray && !is_pad)
      is_rgba = 1;
    if ((!strstr(map,"c") || !strstr(map,"m") || !strstr(map,"y") || !strstr(map,"k")) && !is_rgba && !is_gray && !is_pad)
      is_cmyka = 1;
    if (!strstr(map,"i") && !is_rgba && !is_cmyka && !is_pad)
      is_gray = 1;
    if (!strstr(map,"p") && !is_rgba && !is_cmyka && !is_gray)
      is_pad = 1;
    if (is_cmyka + is_rgba + is_gray + is_pad != 1)
    {
      fprintf(stderr, THIS_SOLVER ": List 'image' has improper 'pixel_map'\n");
      rerror(THIS_SOLVER ": Cannot continue!\n");
    }

    // fix the size
    if (height == -1)
      height = lheight;
    if (width == -1)
      width = lwidth;

    if (qrange==255)
    {
      pixels_8 = (unsigned char *) GC_malloc(height * width * npf * sizeof(unsigned char));
      for (i=0; i<height; i++)
      {
        for (j=0; j<width; j++)
        {
          for (k=0; k<npf; k++)
          {
//             fprintf(stderr, "qrange = %g", qrange);

            if (pixel_storage == 'd')
            {
//               fprintf(stderr, "Mdr0( pixel_data[k],i,j) = %g", Mdr0( pixel_data[k],i,j));
//               fprintf(stderr, " -> %i\n", (unsigned char)(qrange * Mdr0( pixel_data[k],i,j)));
              pixels_8[i*npf*width + npf*j + k] = (unsigned char)(qrange * Mdr0( pixel_data[k],i,j));
            }
            else if (pixel_storage == 'i')
            {
//               fprintf(stderr, "Mdi0( pixel_data[k],i,j) = %i", Mdi0( pixel_data[k],i,j));
//               fprintf(stderr, " -> %i\n", (unsigned char)Mdi0( pixel_data[k],i,j));
              pixels_8[i*npf*width + npf*j + k] = (unsigned char) Mdi0(pixel_data[k],i,j);
            }
          }
        }
      }
      int idx = (int) MagickGetIteratorIndex(magick_wand[default_wand_idx]);
      if (idx != -1)
      {
//         fprintf(stderr, "MagickImportImagePixels: w=%i,h=%i, x=%i,y=%i\n", width, height,x,y);
        // image already exists in the default wand
        MagickImportImagePixels(magick_wand[default_wand_idx], (size_t) x, (size_t) y,
                                (size_t) width, (size_t) height, map, CharPixel, (void *)pixels_8);
      }
      else
      {
        if (x==0 && y==0)
        {
//           fprintf(stderr, "MagickConstituteImage\n");
          // we need to create an image using the provided info
          MagickConstituteImage(magick_wand[default_wand_idx], (size_t) width, (size_t) height,
                                map, CharPixel, (void *)pixels_8);
        }
        else
        {
//           fprintf(stderr, "MagickImportImagePixels-2\n");
          PixelWand *bkgrnd = NewPixelWand();
          PixelSetRed(bkgrnd, qrange);
          PixelSetGreen(bkgrnd, qrange);
          PixelSetBlue(bkgrnd, qrange);
          PixelSetOpacity(bkgrnd, qrange);

          // create new image on new wand
          MagickNewImage(magick_wand[default_wand_idx], (size_t) (y+width), (size_t) (x+height),bkgrnd);

          bkgrnd = DestroyPixelWand(bkgrnd);

          // image already exists in the default wand
          MagickImportImagePixels(magick_wand[default_wand_idx], (size_t) y, (size_t) x,
                                  (size_t) width, (size_t) height, map, CharPixel, (void *)pixels_8);
        }
      }

      GC_free(pixels_8);
    }
    else if (qrange==65535)
    {
      pixels_16 = (short unsigned int *) GC_malloc(height * width * npf * sizeof(short unsigned int));

      for (j=0; j<height; j++)
      {
        for (i=0; i<width; i++)
        {
          for (k=0; k<npf; k++)
          {
            if (pixel_storage == 'd')
              pixels_16[i*npf*width + npf*j + k] = (unsigned short int)(qrange * Mdr0(pixel_data[k],i,j));
            else if (pixel_storage == 'i')
              pixels_16[i*npf*width + npf*j + k] = (unsigned short int) Mdi0(pixel_data[k],i,j) > (unsigned short int) qrange
              ? (unsigned short int) qrange : (unsigned short int) Mdi0(pixel_data[k],i,j);
          }
        }
      }
      MagickImportImagePixels(magick_wand[default_wand_idx], (size_t) x, (size_t) y,
                              (size_t) width, (size_t) height, map, ShortPixel, (void *)pixels_16);

      GC_free(pixels_16);
    }

    ent_Clean (e1);
    ent_Clean (en);

    rent = ent_Create ();
    ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
    ent_type (rent) = MATRIX_DENSE_REAL;
    return (rent);
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //
  // extract pixel arrays from top image in default magick wand
  //
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": First argument '[y1,y2]' must be 2-vector for column range!\n");
  MDR *w = ent_data(e1);
  y = mdiV0(w,0) - 1;
  width = mdiV0(w,SIZE(w)-1) - y;

  // second argument [x,h] or [x:h]
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Second argument '[x1,x2]' must be 2-vector for row range!\n");
  MDR *h = ent_data(e2);
  x = mdiV0(h,0) - 1;
  height = mdiV0(h,SIZE(h)-1) - x;

  // third argument: pixel map
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_STRING)
    rerror(THIS_SOLVER ": Third argument 'pixel_map' must be string.\n");
  map = class_char_pointer (e3);
  if (!map)
  {
    fprintf(stderr, THIS_SOLVER ": Third argument 'pixel_map' must be string.\n");
    rerror(THIS_SOLVER ": Cannot continue!\n");
  }
  if (strlen(map)<1 || strlen(map)>4)
  {
    fprintf(stderr, THIS_SOLVER ": Third argument 'pixel_map' can be between 1 and 4 character long!\n");
    rerror(THIS_SOLVER ": Cannot continue!\n");
  }
  npf = strlen(map);

  // quantum range
  if (nargs == 4)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == MATRIX_DENSE_REAL)
    {
      ddummy = class_double (e4);
      if (ddummy==8 || ddummy==255 || ddummy==1)
        qrange = 255;
      else if (ddummy==16 || ddummy==65535 || ddummy==2)
        qrange = 65535;
    }
  }

  // create necessary arrays
  for (k=0; k<npf; k++)
    pixel_data[k] = mdi_Create(height,width);

  if (qrange==255)
  {
    pixels_8 = (unsigned char *) GC_malloc(height * width * npf * sizeof(unsigned char));

    MagickExportImagePixels(magick_wand[default_wand_idx], (size_t) y, (size_t) x,
                            (size_t) width, (size_t) height, map, CharPixel, (void *)pixels_8);
    for (i=0; i<height; i++)
    {
      for (j=0; j<width; j++)
      {
        for (k=0; k<npf; k++)
        {
          Mdi0(pixel_data[k],i,j) = (int) pixels_8[i*npf*width + npf*j + k];
        }
      }
    }

    GC_free(pixels_8);
  }
  else if (qrange==65535)
  {
    pixels_16 = (unsigned short int *) GC_malloc(height * width * npf * sizeof(unsigned short int));
    MagickExportImagePixels(magick_wand[default_wand_idx], (size_t) y, (size_t) x,
                            (size_t) width, (size_t) height, map, ShortPixel, (void *)pixels_16);
    for (i=0; i<height; i++)
    {
      for (j=0; j<width; j++)
      {
        for (k=0; k<npf; k++)
        {
          Mdi0(pixel_data[k],i,j) = (int) pixels_16[i*npf*width + npf*j + k];
        }
      }
    }

    GC_free(pixels_16);
  }

  //
  // return the result as a list
  //
  Btree *bw = btree_Create ();
  // x:
  Ent *R = ent_Create ();
  ent_data (R) = mdi_CreateScalar(x+1);
  ent_type (R) = MATRIX_DENSE_REAL;
  install (bw, ("x"), R);

  // y:
  Ent *S = ent_Create ();
  ent_data (S) = mdi_CreateScalar(y+1);
  ent_type (S) = MATRIX_DENSE_REAL;
  install (bw, ("y"), S);

  // height:
  Ent *T = ent_Create ();
  ent_data (T) = mdi_CreateScalar(height);
  ent_type (T) = MATRIX_DENSE_REAL;
  install (bw, ("height"), T);

  // width:
  Ent *U = ent_Create ();
  ent_data (U) = mdi_CreateScalar(width);
  ent_type (U) = MATRIX_DENSE_REAL;
  install (bw, ("width"), U);

  // pixel_map:
  Ent *V = ent_Create ();
  ent_data (V) = mds_CreateScalar(map);
  ent_type (V) = MATRIX_DENSE_STRING;
  install (bw, ("pixel_map"), V);

  // sublist: pixel
  Btree *bw2 = btree_Create ();
  for (k=0;k<npf;k++)
  {
    // map contains initials of channels. we expect these to appear in the list pixel
    channel[0] = map[k];
    channel[1] = '\0';

    Ent *W = ent_Create ();
    ent_data (W) = pixel_data[k];
    ent_type (W) = MATRIX_DENSE_REAL;
    install (bw2, (channel), W);
  }
  Ent *Z = ent_Create ();
  ent_data (Z) = bw2;
  ent_type (Z) = BTREE;
  install (bw, ("pixel"), Z);

  // qrange:
  Ent *Q = ent_Create ();
  ent_data (Q) = mdi_CreateScalar(qrange);
  ent_type (Q) = MATRIX_DENSE_REAL;
  install (bw, ("qrange"), Q);

  // last minute cleanup
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  rent = ent_Create ();
  ent_data (rent) = bw;
  ent_type (rent) = BTREE;
  return rent;
}


//
// evaluate operator on image
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".f"
Ent *
ent_im_wand_f (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *es=0, *rent;
  char *c=0, *expr=0;
  MDR *params=0;

  int ar=0, ac=0, ag=0, am=0, ab=0, ay=0, ao=0, aB=0, aM=0, is_eval=0, is_func=0, is_fx=0, len=0;

  MagickEvaluateOperator op=UndefinedEvaluateOperator;
  MagickFunction function=UndefinedFunction;

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs == 0)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Evaluate operator/function/expression on active image from default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(operator[,parameters][,channel]),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'operator' is string containing the operation,\n");
    fprintf (stdout,
             THIS_SOLVER ":   'value' is real scalar parameter of the operator, while\n");
    fprintf (stdout,
             THIS_SOLVER ":   'channel' is the color channel to which the operator is applied.\n");
    rerror ( THIS_SOLVER ": requires 1,2 or 3 arguments");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    rerror(THIS_SOLVER ": First argument 'operator' must be string!");
  expr = class_char_pointer(e1);

  if (!strcmp(expr,"+"))
  {
    is_eval = 1;
    op = AddEvaluateOperator;
  }
  else if (!strcmp(expr,"and"))
  {
    is_eval = 1;
    op = AndEvaluateOperator;
  }
  else if (!strcmp(expr,"/") || !strcmp(expr,"divide"))
  {
    is_eval = 1;
    op = DivideEvaluateOperator;
  }
  else if (!strcmp(expr,"<<") || !strcmp(expr,"shiftleft"))
  {
    is_eval = 1;
    op = LeftShiftEvaluateOperator;
  }
  else if (!strcmp(expr,"max"))
  {
    is_eval = 1;
    op = MaxEvaluateOperator;
  }
  else if (!strcmp(expr,"min"))
  {
    is_eval = 1;
    op = MinEvaluateOperator;
  }
  else if (!strcmp(expr,"*") || !strcmp(expr,"multiply"))
  {
    is_eval = 1;
    op = MultiplyEvaluateOperator;
  }
  else if (!strcmp(expr,"or"))
  {
    is_eval = 1;
    op = OrEvaluateOperator;
  }
  else if (!strcmp(expr,">>") || !strcmp(expr,"shiftright"))
  {
    is_eval = 1;
    op = RightShiftEvaluateOperator;
  }
  else if (!strcmp(expr,"set"))
  {
    is_eval = 1;
    op = SetEvaluateOperator;
  }
  else if (!strcmp(expr,"-") || !strcmp(expr,"subtract") || !strcmp(expr,"minus"))
  {
    is_eval = 1;
    op = SubtractEvaluateOperator;
  }
  else if (!strcmp(expr,"xor"))
  {
    is_eval = 1;
    op = XorEvaluateOperator;
  }
  else if (!strcmp(expr,"pow"))
  {
    is_eval = 1;
    op = PowEvaluateOperator;
  }
  else if (!strcmp(expr,"log"))
  {
    is_eval = 1;
    op = LogEvaluateOperator;
  }
  else if (!strcmp(expr,"threshold"))
  {
    is_eval = 1;
    op = ThresholdEvaluateOperator;
  }
  else if (!strcmp(expr,"thresholdblack"))
  {
    is_eval = 1;
    op = ThresholdBlackEvaluateOperator;
  }
  else if (!strcmp(expr,"thresholdwhite"))
  {
    is_eval = 1;
    op = ThresholdWhiteEvaluateOperator;
  }
  else if (!strcmp(expr,"gaussian") || !strcmp(expr,"noisegaussian"))
  {
    is_eval = 1;
    op = GaussianNoiseEvaluateOperator;
  }
  else if (!strcmp(expr,"impulse") || !strcmp(expr,"noiseimpulse"))
  {
    is_eval = 1;
    op = ImpulseNoiseEvaluateOperator;
  }
  else if (!strcmp(expr,"laplace") || !strcmp(expr,"noiselaplace"))
  {
    is_eval = 1;
    op = LaplacianNoiseEvaluateOperator;
  }
  else if (!strcmp(expr,"multiplicative") || !strcmp(expr,"noisemultiplicative"))
  {
    is_eval = 1;
    op = MultiplicativeNoiseEvaluateOperator;
  }
  else if (!strcmp(expr,"poisson") || !strcmp(expr,"noisepoisson"))
  {
    is_eval = 1;
    op = PoissonNoiseEvaluateOperator;
  }
  else if (!strcmp(expr,"uniform") || !strcmp(expr,"noiseuniform"))
  {
    is_eval = 1;
    op = UniformNoiseEvaluateOperator;
  }
  else if (!strcmp(expr,"cos"))
  {
    is_eval = 1;
    op = CosineEvaluateOperator;
  }
  else if (!strcmp(expr,"sin"))
  {
    is_eval = 1;
    op = SineEvaluateOperator;
  }
  else if (!strcmp(expr,"+mod") || !strcmp(expr,"addmodulus"))
  {
    is_eval = 1;
    op = AddModulusEvaluateOperator;
  }
  else if (!strcmp(expr,"mean"))
  {
    is_eval = 1;
    op = MeanEvaluateOperator;
  }
  else if (!strcmp(expr,"abs"))
  {
    is_eval = 1;
    op = AbsEvaluateOperator;
  }
  else if (!strcmp(expr,"exp"))
  {
    is_eval = 1;
    op = ExponentialEvaluateOperator;
  }
  else if (!strcmp(expr,"median"))
  {
    is_eval = 1;
    op = MedianEvaluateOperator;
  }
  //     else if (!strcmp(c,"sum"))
    //       op = SumEvaluateOperator;
  else if (!strcmp(c,"poly"))
  {
    is_func = 1;
    function = PolynomialFunction;
  }
//   else if (!strcmp(c,"sinf"))
//   {
//     is_func = 1;
//     function = SinusoidFunction;
//   }
  else if (!strcmp(c,"asin"))
  {
    is_func = 1;
    function = ArcsinFunction;
  }
  else if (!strcmp(c,"atan"))
  {
    is_func = 1;
    function = ArctanFunction;
  }
  else
  {
    // assume it is expression and just pass it along
    is_fx = 1;
  }


  //
  // process channel info if provided
  //
  if (nargs > 1)
  {
    es = bltin_get_ent (args[nargs-1]);
    if (ent_type(es) == MATRIX_DENSE_STRING)
    {
      c = class_char_pointer(es);
      if (c)
      {
        //
        // process channels if necessary
        //
        ar = (strstr(c, "r") != 0);  // RedChannel
        ac = (strstr(c, "c") != 0);  // CyanChannel
        ag = (strstr(c, "g") != 0);  // GreenChannel
        am = (strstr(c, "m") != 0);  // MagentaChannel
        ab = (strstr(c, "b") != 0);  // BlueChannel
        ay = (strstr(c, "y") != 0);  // YellowChannel
        ao = (strstr(c, "o") != 0 || strstr(c, "a") != 0);  // OpacityChannel
        aB = (strstr(c, "B") != 0 || strstr(c, "k") != 0);  // BlackChannel
        aM = (strstr(c, "M") != 0 || strstr(c, "t") != 0);  // MatteChannel
      }
    }
    ent_Clean(es);
  }

  if (nargs >= 2 && !is_fx)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) == MATRIX_DENSE_REAL)
    {
      params = ent_data(e2);
      if (params->type != RLAB_TYPE_DOUBLE)
        rerror(THIS_SOLVER ": Second argument, if numeric, must be real vector!");
      len = params->nrow * params->ncol;
    }
  }

  if (is_eval)
  {
    //
    // call evaluate
    //
    if (ar)
      MagickEvaluateImageChannel(magick_wand[default_wand_idx], RedChannel, op, MdrV0(params,0));
    else if (ac)
      MagickEvaluateImageChannel(magick_wand[default_wand_idx], CyanChannel, op, MdrV0(params,0));
    else if (ag)
      MagickEvaluateImageChannel(magick_wand[default_wand_idx], GreenChannel, op, MdrV0(params,0));
    else if (am)
      MagickEvaluateImageChannel(magick_wand[default_wand_idx], MagentaChannel, op, MdrV0(params,0));
    else if (ab)
      MagickEvaluateImageChannel(magick_wand[default_wand_idx], BlueChannel, op, MdrV0(params,0));
    else if (ay)
      MagickEvaluateImageChannel(magick_wand[default_wand_idx], YellowChannel, op, MdrV0(params,0));
    else if (ao)
      MagickEvaluateImageChannel(magick_wand[default_wand_idx], OpacityChannel, op, MdrV0(params,0));
    else if (aB)
      MagickEvaluateImageChannel(magick_wand[default_wand_idx], BlackChannel, op, MdrV0(params,0));
    else if (aM)
      MagickEvaluateImageChannel(magick_wand[default_wand_idx], MatteChannel, op, MdrV0(params,0));
    else if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
      MagickEvaluateImage(magick_wand[default_wand_idx], op, MdrV0(params,0));
  }
  else if (is_fx)
  {
    //
    // call fx
    //
    if (ar)
      MagickFxImageChannel(magick_wand[default_wand_idx], RedChannel, expr);
    else if (ac)
      MagickFxImageChannel(magick_wand[default_wand_idx], CyanChannel, expr);
    else if (ag)
      MagickFxImageChannel(magick_wand[default_wand_idx], GreenChannel, expr);
    else if (am)
      MagickFxImageChannel(magick_wand[default_wand_idx], MagentaChannel, expr);
    else if (ab)
      MagickFxImageChannel(magick_wand[default_wand_idx], BlueChannel, expr);
    else if (ay)
      MagickFxImageChannel(magick_wand[default_wand_idx], YellowChannel, expr);
    else if (ao)
      MagickFxImageChannel(magick_wand[default_wand_idx], OpacityChannel, expr);
    else if (aB)
      MagickFxImageChannel(magick_wand[default_wand_idx], BlackChannel, expr);
    else if (aM)
      MagickFxImageChannel(magick_wand[default_wand_idx], MatteChannel, expr);
    else if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
      MagickFxImage(magick_wand[default_wand_idx], expr);
  }
  else if (is_func)
  {
    if (ar)
      MagickFunctionImageChannel(magick_wand[default_wand_idx], RedChannel, function, (size_t) len, MDRPTR(params));
    else if (ac)
      MagickFunctionImageChannel(magick_wand[default_wand_idx], CyanChannel, function, (size_t) len, MDRPTR(params));
    else if (ag)
      MagickFunctionImageChannel(magick_wand[default_wand_idx], GreenChannel, function, (size_t) len, MDRPTR(params));
    else if (am)
      MagickFunctionImageChannel(magick_wand[default_wand_idx], MagentaChannel, function, (size_t) len, MDRPTR(params));
    else if (ab)
      MagickFunctionImageChannel(magick_wand[default_wand_idx], BlueChannel, function, (size_t) len, MDRPTR(params));
    else if (ay)
      MagickFunctionImageChannel(magick_wand[default_wand_idx], YellowChannel, function, (size_t) len, MDRPTR(params));
    else if (ao)
      MagickFunctionImageChannel(magick_wand[default_wand_idx], OpacityChannel, function, (size_t) len, MDRPTR(params));
    else if (aB)
      MagickFunctionImageChannel(magick_wand[default_wand_idx], BlackChannel, function, (size_t) len, MDRPTR(params));
    else if (aM)
      MagickFunctionImageChannel(magick_wand[default_wand_idx], MatteChannel, function, (size_t) len, MDRPTR(params));
    else if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
      MagickFunctionImage(magick_wand[default_wand_idx], function, (size_t) len, MDRPTR(params));
  }

  ent_Clean(e1);
  ent_Clean(e2);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}


#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".specs"
Ent *
ent_im_wand_image_specs(int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  char *spec=0;
  int i;
  double d;
  MDR *w=0;

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs != 0 && nargs !=2

  )
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Specify or retrive specifications of active image in default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ": (1) specs  = " THIS_SOLVER "(),\n");
    fprintf (stdout,
             THIS_SOLVER ": (2) " THIS_SOLVER "(spec,val),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'spec' is string, e.g., from list 'specs',\n");
    fprintf (stdout,
             THIS_SOLVER ":   'val' is the value of that perticular specification in the active image.,\n");
    fprintf (stdout,
             THIS_SOLVER ": See manual for more info.\n");
    rerror ( THIS_SOLVER ": requires none or 2 arguments");
  }

  if (nargs == 2)
  {
    //
    // what specification are we talking about
    //
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) != MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": First argument 'spec' must be string!");
    spec = class_char_pointer(e1);
    if (!spec)
      rerror(THIS_SOLVER ": First argument 'spec' must be string!");

    rent = ent_Create();
    //
    // process the value according to the specification
    //
    if (!strcmp(spec, "alpha"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be integer or real scalar!");
      d = class_double(e2);
      if (d<0 || d>1)
        rerror(THIS_SOLVER ": Second argument 'val' must be integer or real scalar!");
      if (d==1)
        ent_data(rent) = mdi_CreateScalar( MagickSetImageAlphaChannel(magick_wand[default_wand_idx],
                                                                     ActivateAlphaChannel) );
      else if (d==0)
        ent_data(rent) = mdi_CreateScalar( MagickSetImageAlphaChannel(magick_wand[default_wand_idx],
                                                                      DeactivateAlphaChannel) );
      else
        ent_data(rent) = mdi_CreateScalar( MagickSetImageOpacity(magick_wand[default_wand_idx], d) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "delay"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be integer scalar!");
      i = (int) class_double(e2);

      ent_data(rent) = mdi_CreateScalar( MagickSetImageDelay(magick_wand[default_wand_idx],
                                                             (size_t) i) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "fuzz"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real scalar!");
      d = (int) class_double(e2);

      ent_data(rent) = mdi_CreateScalar( MagickSetImageFuzz(magick_wand[default_wand_idx], d) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "depth"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be integer scalar!");
      i = (int) class_double(e2);

      ent_data(rent) = mdi_CreateScalar( MagickSetImageDepth(magick_wand[default_wand_idx],
                                                             (size_t) i) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "gamma"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real scalar!");
      d = (int) class_double(e2);

      ent_data(rent) = mdi_CreateScalar( MagickSetImageGamma(magick_wand[default_wand_idx], d) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "page"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be integer vector!");
      w = class_matrix_real(e2);
      if (!w)
        rerror(THIS_SOLVER ": Second argument 'val' must be integer vector!");
      if (w->nrow * w->ncol != 4)
        rerror(THIS_SOLVER ": Second argument 'val' must be integer vector!");

      size_t width  = (size_t) mdrV0(w,0);
      size_t height = (size_t) mdrV0(w,1);
      ssize_t x = (ssize_t) mdrV0(w,2);
      ssize_t y = (ssize_t) mdrV0(w,3);

      ent_data(rent) = mdi_CreateScalar( MagickSetImagePage(magick_wand[default_wand_idx],
                                                            width, height, x, y) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "background"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
      if (!w)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      if (w->nrow * w->ncol != 3 && w->nrow * w->ncol != 4 )
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");

      PixelWand *pw = NewPixelWand();
      PixelSetRed     (pw, (double) MdrV0(w,0));
      PixelSetGreen   (pw, (double) MdrV0(w,1));
      PixelSetBlue    (pw, (double) MdrV0(w,2));
      if (w->nrow * w->ncol == 4)
        PixelSetOpacity (pw, (double) MdrV0(w,3));
      else
        PixelSetOpacity (pw, (double) 0);
      ent_data(rent) = mdi_CreateScalar( MagickSetImageBackgroundColor(magick_wand[default_wand_idx],
                                                                       pw) );
      pw = DestroyPixelWand(pw);
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "bias"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real scalar!");
      d = (int) class_double(e2);

      ent_data(rent) = mdi_CreateScalar( MagickSetImageBias(magick_wand[default_wand_idx], d) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "blue_primary"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
      if (!w)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      if (w->nrow * w->ncol != 2)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");

      ent_data(rent) = mdi_CreateScalar( MagickSetImageBluePrimary(magick_wand[default_wand_idx],
                                                                   mdrV0(w,0), mdrV0(w,1)) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "red_primary"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
      if (!w)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      if (w->nrow * w->ncol != 2)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");

      ent_data(rent) = mdi_CreateScalar( MagickSetImageRedPrimary(magick_wand[default_wand_idx],
                                                                  mdrV0(w,0), mdrV0(w,1)) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "green_primary"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
      if (!w)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      if (w->nrow * w->ncol != 2)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");

      ent_data(rent) = mdi_CreateScalar( MagickSetImageGreenPrimary(magick_wand[default_wand_idx],
                                                                    mdrV0(w,0), mdrV0(w,1)) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "border_color"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
      if (!w)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      if (w->nrow * w->ncol != 3 && w->nrow * w->ncol != 4 )
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");

      PixelWand *pw = NewPixelWand();
      PixelSetRed     (pw, (double) MdrV0(w,0));
      PixelSetGreen   (pw, (double) MdrV0(w,1));
      PixelSetBlue    (pw, (double) MdrV0(w,2));
      if (w->nrow * w->ncol == 4)
        PixelSetOpacity (pw, (double) MdrV0(w,3));
      else
        PixelSetOpacity (pw, (double) 0);
      ent_data(rent) = mdi_CreateScalar( MagickSetImageBorderColor(magick_wand[default_wand_idx], pw) );
      pw = DestroyPixelWand(pw);
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "color"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
      if (!w)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      if (w->nrow * w->ncol != 3 && w->nrow * w->ncol != 4 )
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");

      PixelWand *pw = NewPixelWand();
      PixelSetRed     (pw, (double) MdrV0(w,0));
      PixelSetGreen   (pw, (double) MdrV0(w,1));
      PixelSetBlue    (pw, (double) MdrV0(w,2));
      if (w->nrow * w->ncol == 4)
        PixelSetOpacity (pw, (double) MdrV0(w,3));
      else
        PixelSetOpacity (pw, (double) 0);
      ent_data(rent) = mdi_CreateScalar( MagickSetImageColor(magick_wand[default_wand_idx], pw) );
      pw = DestroyPixelWand(pw);
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "color_space"))
    {
      ColorspaceType cs;
      char *val=0;

      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");
      val = class_char_pointer(e2);
      if (!val)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");

      if (!strcmp(val,"rgb" ))
        cs = RGBColorspace;
      else if (!strcmp(val,"gray"))
        cs = GRAYColorspace;
      else if (!strcmp(val,"transparent"))
        cs = TransparentColorspace;
      else if (!strcmp(val,"ohta"))
        cs = OHTAColorspace;
      else if (!strcmp(val,"lab"))
        cs = LabColorspace;
      else if (!strcmp(val,"xyz"))
        cs = XYZColorspace;
      else if (!strcmp(val,"ycbcr"))
        cs = YCbCrColorspace;
      else if (!strcmp(val,"ycc"))
        cs = YCCColorspace;
      else if (!strcmp(val,"yiq"))
        cs = YIQColorspace;
      else if (!strcmp(val,"ypbpr"))
        cs = YPbPrColorspace;
      else if (!strcmp(val,"yuv"))
        cs = YUVColorspace;
      else if (!strcmp(val,"cmyk"))
        cs = CMYKColorspace;
      else if (!strcmp(val,"srgb"))
        cs = sRGBColorspace;
      else if (!strcmp(val,"hsb"))
        cs = HSBColorspace;
      else if (!strcmp(val,"hsl"))
        cs = HSLColorspace;
      else if (!strcmp(val,"hwb"))
        cs = HWBColorspace;
      else if (!strcmp(val,"rec601luma"))
        cs = Rec601LumaColorspace;
      else if (!strcmp(val,"rec601ycbcr"))
        cs = Rec601YCbCrColorspace;
      else if (!strcmp(val,"rec709luma"))
        cs = Rec709LumaColorspace;
      else if (!strcmp(val,"rec709ycbcr"))
        cs = Rec709YCbCrColorspace;
      else if (!strcmp(val,"log"))
        cs = LogColorspace;
      else if (!strcmp(val,"cmy"))
        cs = CMYColorspace;
      else
        cs = RGBColorspace;

      ent_data(rent) = mdi_CreateScalar( MagickSetImageColorspace( magick_wand[default_wand_idx], cs) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "filename"))
    {
      char *val=0;

      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");
      val = class_char_pointer(e2);
      if (!val)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");

      ent_data(rent) = mdi_CreateScalar( MagickSetImageFilename( magick_wand[default_wand_idx], val) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "format"))
    {
      char *val=0;

      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");
      val = class_char_pointer(e2);
      if (!val)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");

      ent_data(rent) = mdi_CreateScalar( MagickSetImageFormat( magick_wand[default_wand_idx], val) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "gravity"))
    {
      // Oh gravity! It's pulling me!
      char *val=0;
      GravityType gravity;

      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");
      val = class_char_pointer(e2);
      if (!val)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");

      if (!strcmp(val,"nw"))
        gravity=NorthWestGravity;
      else if (!strcmp(val,"n"))
        gravity=NorthGravity;
      else if (!strcmp(val,"ne"))
        gravity=NorthEastGravity;
      else if (!strcmp(val,"w"))
        gravity=WestGravity;
      else if (!strcmp(val,"center"))
        gravity=CenterGravity;
      else if (!strcmp(val,"e"))
        gravity=EastGravity;
      else if (!strcmp(val,"sw"))
        gravity=SouthWestGravity;
      else if (!strcmp(val,"s"))
        gravity=SouthGravity;
      else if (!strcmp(val,"se"))
        gravity=SouthEastGravity;
      else
        gravity=ForgetGravity;

      ent_data(rent) = mdi_CreateScalar( MagickSetImageGravity( magick_wand[default_wand_idx], gravity) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "interlace"))
    {
      // Oh gravity! It's pulling me!
      char *val=0;
      InterlaceType interlace;

      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");
      val = class_char_pointer(e2);
      if (!val)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");

      if (!strcmp(val,"line"))
        interlace=LineInterlace;
      else if (!strcmp(val,"plane"))
        interlace=PlaneInterlace;
      else if (!strcmp(val,"part"))
        interlace=PartitionInterlace;
      else if (!strcmp(val,"gif"))
        interlace=GIFInterlace;
      else if (!strcmp(val,"jpeg"))
        interlace=JPEGInterlace;
      else if (!strcmp(val,"png"))
        interlace=PNGInterlace;
      else
        interlace=NoInterlace;

      ent_data(rent) = mdi_CreateScalar( MagickSetImageInterlaceScheme( magick_wand[default_wand_idx],
                                                                        interlace ) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "interpolate"))
    {
      // Oh gravity! It's pulling me!
      char *val=0;
      InterpolatePixelMethod interpolate;

      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");
      val = class_char_pointer(e2);
      if (!val)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");

      if (!strcmp(val,"bicubic"))
        interpolate=BicubicInterpolatePixel;
      else if (!strcmp(val,"bilinear"))
        interpolate=BilinearInterpolatePixel;
      else if (!strcmp(val,"filter"))
        interpolate=FilterInterpolatePixel;
      else if (!strcmp(val,"int"))
        interpolate=IntegerInterpolatePixel;
      else if (!strcmp(val,"nn"))
        interpolate=NearestNeighborInterpolatePixel;
      else if (!strcmp(val,"spline"))
        interpolate=SplineInterpolatePixel;
      else
        interpolate=AverageInterpolatePixel;

      ent_data(rent) = mdi_CreateScalar( MagickSetImageInterpolateMethod( magick_wand[default_wand_idx],
                                                                          interpolate ) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "matte"))
    {
      MagickBooleanType matte = MagickFalse;

      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be integer or real scalar!");
      w = class_matrix_real (e2);
      if (w->nrow * w->ncol != 1 && w->nrow * w->ncol != 3 && w->nrow * w->ncol != 4)
        rerror(THIS_SOLVER ": Second argument 'val' must be integer or real scalar!");

      if (w->nrow * w->ncol == 1)
      {
        if (mdrV0(w,0))
          matte = MagickTrue;

        ent_data(rent) = mdi_CreateScalar( MagickSetImageMatte(magick_wand[default_wand_idx],matte) );
      }
      else
      {
        PixelWand *pw = NewPixelWand();
        PixelSetRed     (pw, (double) MdrV0(w,0));
        PixelSetGreen   (pw, (double) MdrV0(w,1));
        PixelSetBlue    (pw, (double) MdrV0(w,2));
        if (w->nrow * w->ncol == 4)
          PixelSetOpacity (pw, (double) MdrV0(w,3));
        else
          PixelSetOpacity (pw, (double) 0);

        ent_data(rent) = mdi_CreateScalar( MagickSetImageMatteColor(magick_wand[default_wand_idx], pw) );

        pw = DestroyPixelWand(pw);
        ent_type(rent) = MATRIX_DENSE_REAL;
      }

      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "resolution"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
      if (!w)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      if (w->nrow * w->ncol != 2)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");

      ent_data(rent) = mdi_CreateScalar( MagickSetImageResolution(magick_wand[default_wand_idx],
                                                                  mdrV0(w,0), mdrV0(w,1)) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "type"))
    {
      // Oh gravity! It's pulling me!
      char *val=0;
      ImageType imagetype=UndefinedType;

      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");
      val = class_char_pointer(e2);
      if (!val)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");

      if (!strcmp(val,"bilevel"))
        imagetype=BilevelType;
      else if (!strcmp(val,"grayscale"))
        imagetype=GrayscaleType;
      else if (!strcmp(val,"grayscalematte"))
        imagetype=GrayscaleMatteType;
      else if (!strcmp(val,"palette"))
        imagetype=PaletteType;
      else if (!strcmp(val,"palettematte"))
        imagetype=PaletteMatteType;
      else if (!strcmp(val,"truecolor"))
        imagetype=TrueColorType;
      else if (!strcmp(val,"truecolormatte"))
        imagetype=TrueColorMatteType;
      else if (!strcmp(val,"colorseparation"))
        imagetype=ColorSeparationType;
      else if (!strcmp(val,"colorseparationmatte"))
        imagetype=ColorSeparationMatteType;
      else if (!strcmp(val,"optimize"))
        imagetype=OptimizeType;
      else if (!strcmp(val,"palettebilevelmatte"))
        imagetype=PaletteBilevelMatteType;
      else
        imagetype=UndefinedType;

      ent_data(rent) = mdi_CreateScalar( MagickSetImageType( magick_wand[default_wand_idx],
                                                                          imagetype ) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "resolution_unit"))
    {
      // Give me metric system or give me death!
      char *val=0;
      ResolutionType units=UndefinedResolution;

      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");
      val = class_char_pointer(e2);
      if (!val)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");

      if (!strcmp(val,"ppi"))
        units=PixelsPerInchResolution;
      else if (!strcmp(val,"ppcm"))
        units=PixelsPerCentimeterResolution;
      else
        units=UndefinedResolution;

      ent_data(rent) = mdi_CreateScalar( MagickSetImageUnits( magick_wand[default_wand_idx],
                                                              units ) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "white_point"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
      if (!w)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      if (w->nrow * w->ncol != 2)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");

      ent_data(rent) = mdi_CreateScalar( MagickSetImageWhitePoint(magick_wand[default_wand_idx],
                                                                  mdrV0(w,0), mdrV0(w,1)) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else
    {
      fprintf(stderr, THIS_SOLVER  ": Option not recognized. Consult manual!\n");
      ent_data(rent) = mdi_Create(0,0);
      ent_type(rent) = MATRIX_DENSE_REAL;
    }

    ent_Clean (e1);
    ent_Clean (e2);
    return (rent);
  }

  //
  // no parameters: go over all specifications that exist
  //

  //
  // return the result as a list
  //
  Btree *bw = btree_Create ();

  Ent *E;

  // alpha:
  E = ent_Create ();
  ent_data (E) = mdi_CreateScalar( (int) MagickGetImageAlphaChannel(magick_wand[default_wand_idx]) );
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("alpha"), E);

  // red primary:
  MDR *r = mdr_Create(1,2);
  MagickGetImageRedPrimary(magick_wand[default_wand_idx], &MdrV0(r,0), &MdrV0(r,1));
  E = ent_Create ();
  ent_data (E) = r;
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("red_primary"), E);

  // green primary:
  MDR *g = mdr_Create(1,2);
  MagickGetImageGreenPrimary(magick_wand[default_wand_idx], &MdrV0(g,0), &MdrV0(g,1));
  E = ent_Create ();
  ent_data (E) = g;
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("green_primary"), E);

  // blue primary:
  MDR *b = mdr_Create(1,2);
  MagickGetImageBluePrimary(magick_wand[default_wand_idx], &MdrV0(b,0), &MdrV0(b,1));
  E = ent_Create ();
  ent_data (E) = b;
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("blue_primary"), E);

  // colorspace
  ColorspaceType cs = MagickGetImageColorspace( magick_wand[default_wand_idx] );
  E = ent_Create ();
  MDS *sp=0;
  switch (cs)
  {
    case RGBColorspace:
      sp = mds_CreateScalar("rgb");
      break;

    case GRAYColorspace:
      sp = mds_CreateScalar("gray");
      break;

    case TransparentColorspace:
      sp = mds_CreateScalar("transparent");
      break;

    case OHTAColorspace:
      sp = mds_CreateScalar("ohta");
      break;

    case LabColorspace:
      sp = mds_CreateScalar("lab");
      break;

    case XYZColorspace:
      sp = mds_CreateScalar("xyz");
      break;

    case YCbCrColorspace:
      sp = mds_CreateScalar("ycbcr");
      break;

    case YCCColorspace:
      sp = mds_CreateScalar("ycc");
      break;

    case YIQColorspace:
      sp = mds_CreateScalar("yiq");
      break;

    case YPbPrColorspace:
      sp = mds_CreateScalar("ypbpr");
      break;

    case YUVColorspace:
      sp = mds_CreateScalar("yuv");
      break;

    case CMYKColorspace:
      sp = mds_CreateScalar("cmyk");
      break;

    case sRGBColorspace:
      sp = mds_CreateScalar("srgb");
      break;

    case HSBColorspace:
      sp = mds_CreateScalar("hsb");
      break;

    case HSLColorspace:
      sp = mds_CreateScalar("hsl");
      break;

    case HWBColorspace:
      sp = mds_CreateScalar("hwb");
      break;

    case Rec601LumaColorspace:
      sp = mds_CreateScalar("rec601luma");
      break;

    case Rec601YCbCrColorspace:
      sp = mds_CreateScalar("rec601ycbcr");
      break;

    case Rec709LumaColorspace:
      sp = mds_CreateScalar("rec709luma");
      break;

    case Rec709YCbCrColorspace:
      sp = mds_CreateScalar("rec709ycbcr");
      break;

    case LogColorspace:
      sp = mds_CreateScalar("log");
      break;

    case CMYColorspace:
      sp = mds_CreateScalar("cmy");
      break;

    default:
      sp = mds_CreateScalar("undef");
      break;
  }
  ent_data (E) = sp;
  ent_type (E) = MATRIX_DENSE_STRING;
  install (bw, ("color_space"), E);

  // compression
  ColorspaceType ct = MagickGetImageCompression( magick_wand[default_wand_idx] );
  E = ent_Create ();
  MDS *coms=0;
  switch (ct)
  {
    case NoCompression:
      coms = mds_CreateScalar("none");
      break;

    case BZipCompression:
      coms = mds_CreateScalar("bzip");
      break;

    case DXT1Compression:
      coms = mds_CreateScalar("dxt1");
      break;

    case DXT3Compression:
      coms = mds_CreateScalar("dxt3");
      break;

    case DXT5Compression:
      coms = mds_CreateScalar("dxt5");
      break;

    case FaxCompression:
      coms = mds_CreateScalar("fax");
      break;

    case Group4Compression:
      coms = mds_CreateScalar("group4");
      break;

    case JPEGCompression:
      coms = mds_CreateScalar("jpeg");
      break;

    case JPEG2000Compression:
      coms = mds_CreateScalar("jpeg_2000");
      break;

    case LosslessJPEGCompression:
      coms = mds_CreateScalar("jpeg_lossless");
      break;

    case LZWCompression:
      coms = mds_CreateScalar("lzw");
      break;

    case RLECompression:
      coms = mds_CreateScalar("rle");
      break;

    case ZipCompression:
      coms = mds_CreateScalar("zip");
      break;

    case ZipSCompression:
      coms = mds_CreateScalar("zips");
      break;

    case PizCompression:
      coms = mds_CreateScalar("piz");
      break;

    case Pxr24Compression:
      coms = mds_CreateScalar("pxr24");
      break;

    case B44Compression:
      coms = mds_CreateScalar("b44");
      break;

    case B44ACompression:
      coms = mds_CreateScalar("b44a");
      break;

    case LZMACompression:
      coms = mds_CreateScalar("lzma");
      break;

    case JBIG1Compression:
      coms = mds_CreateScalar("jbig1");
      break;

    case JBIG2Compression:
      coms = mds_CreateScalar("jbig2");
      break;

    default:
      coms = mds_CreateScalar("undef");
      break;
  }
  ent_data (E) = coms;
  ent_type (E) = MATRIX_DENSE_STRING;
  install (bw, ("compression"), E);

  // compression quality
  E = ent_Create ();
  ent_data (E) = mdi_CreateScalar( MagickGetImageCompressionQuality(magick_wand[default_wand_idx]) );
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("compression_quality"), E);

  // delay
  E = ent_Create ();
  ent_data (E) = mdi_CreateScalar( MagickGetImageDelay(magick_wand[default_wand_idx]) );
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("delay"), E);

  // depth
  E = ent_Create ();
  ent_data (E) = mdi_CreateScalar( MagickGetImageDepth(magick_wand[default_wand_idx]) );
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("depth"), E);

  // size
  MDR * si = mdi_Create(1,2);
  MdiV0(si,0) = MagickGetImageWidth (magick_wand[default_wand_idx]);
  MdiV0(si,1) = MagickGetImageHeight(magick_wand[default_wand_idx]);
  E = ent_Create ();
  ent_data (E) = si;
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("dim"), E);

//   // endian
//   EndianType endi = MagickGetImageEndian(magick_wand[default_wand_idx]);
//   E = ent_Create ();
//   MDS *end=0;
//   switch (endi)
//   {
//     case LSBEndian:
//       end = mds_CreateScalar("lsb");
//       break;
//
//     case MSBEndian:
//       end = mds_CreateScalar("msb");
//       break;
//
//     default:
//       end = mds_CreateScalar("undef");
//       break;
//   }
//   ent_data (E) = end;
//   ent_type (E) = MATRIX_DENSE_STRING;
//   install (bw, ("endianess"), E);
//   ent_Clean (E);

  // filename
  E = ent_Create();
  ent_data (E) = mds_CreateScalar( MagickGetImageFilename(magick_wand[default_wand_idx]) );
  ent_type (E) = MATRIX_DENSE_STRING;
  install (bw, ("filename"), E);

  // format
  E = ent_Create();
  ent_data (E) = mds_CreateScalar( MagickGetImageFormat(magick_wand[default_wand_idx]) );
  ent_type (E) = MATRIX_DENSE_STRING;
  install (bw, ("format"), E);

  // fuzz
  E = ent_Create();
  ent_data (E) = mdr_CreateScalar( MagickGetImageFuzz(magick_wand[default_wand_idx]) );
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("fuzz"), E);

  // gamma
  E = ent_Create();
  ent_data (E) = mdr_CreateScalar( MagickGetImageGamma(magick_wand[default_wand_idx]) );
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("gamma"), E);

  // gravity
  GravityType gravity = MagickGetImageGravity(magick_wand[default_wand_idx]);
  MDS *grav=0;
  switch (gravity)
  {
    case ForgetGravity:
      grav = mds_CreateScalar("forget");
      break;

    case NorthWestGravity:
      grav = mds_CreateScalar("nw");
      break;

    case NorthGravity:
      grav = mds_CreateScalar("n");
      break;

    case NorthEastGravity:
      grav = mds_CreateScalar("ne");
      break;

    case WestGravity:
      grav = mds_CreateScalar("w");
      break;

    case CenterGravity:
      grav = mds_CreateScalar("center");
      break;

    case EastGravity:
      grav = mds_CreateScalar("e");
      break;

    case SouthWestGravity:
      grav = mds_CreateScalar("sw");
      break;

    case SouthGravity:
      grav = mds_CreateScalar("s");
      break;

    case SouthEastGravity:
      grav = mds_CreateScalar("se");
      break;

    default:
      grav = mds_CreateScalar("undef");
      break;
  }
  E = ent_Create();
  ent_data (E) = grav;
  ent_type (E) = MATRIX_DENSE_STRING;
  install (bw, ("gravity"), E);

  // interlace
  InterlaceType interlace = MagickGetImageInterlaceScheme(magick_wand[default_wand_idx]);
  MDS *intrl=0;
  switch (interlace)
  {
    case NoInterlace:
      intrl = mds_CreateScalar("no");
      break;

    case LineInterlace:
      intrl = mds_CreateScalar("line");
      break;

    case PlaneInterlace:
      intrl = mds_CreateScalar("plane");
      break;

    case PartitionInterlace:
      intrl = mds_CreateScalar("part");
      break;

    case GIFInterlace:
      intrl = mds_CreateScalar("gif");
      break;

    case JPEGInterlace:
      intrl = mds_CreateScalar("jpeg");
      break;

    case PNGInterlace:
      intrl = mds_CreateScalar("png");
      break;

    default:
      intrl = mds_CreateScalar("undef");
      break;
  }
  E = ent_Create();
  ent_data (E) = intrl;
  ent_type (E) = MATRIX_DENSE_STRING;
  install (bw, ("interlace"), E);

  // interpolate
  InterpolatePixelMethod interpolate = MagickGetImageInterpolateMethod(magick_wand[default_wand_idx]);
  MDS *interp=0;
  switch (interpolate)
  {
    case AverageInterpolatePixel:
      interp = mds_CreateScalar("avg");
      break;

    case BicubicInterpolatePixel:
      interp = mds_CreateScalar("bicubic");
      break;

    case BilinearInterpolatePixel:
      interp = mds_CreateScalar("bilinear");
      break;

    case FilterInterpolatePixel:
      interp = mds_CreateScalar("filter");
      break;

    case IntegerInterpolatePixel:
      interp = mds_CreateScalar("int");
      break;

    case NearestNeighborInterpolatePixel:
      interp = mds_CreateScalar("nn");
      break;

    case SplineInterpolatePixel:
      interp = mds_CreateScalar("spline");
      break;

    default:
      interp = mds_CreateScalar("undef");
      break;
  }
  E = ent_Create();
  ent_data (E) = interp;
  ent_type (E) = MATRIX_DENSE_STRING;
  install (bw, ("interpolate"), E);

  // orientation
  InterpolatePixelMethod orientation = MagickGetImageInterpolateMethod(magick_wand[default_wand_idx]);
  MDS *orient=0;
  switch (orientation)
  {
    case TopLeftOrientation:
      orient = mds_CreateScalar("topleft");
      break;

    case TopRightOrientation:
      orient = mds_CreateScalar("topright");
      break;

    case BottomRightOrientation:
      orient = mds_CreateScalar("bottomright");
      break;

    case BottomLeftOrientation:
      orient = mds_CreateScalar("bottomleft");
      break;

    case LeftTopOrientation:
      orient = mds_CreateScalar("lefttop");
      break;

    case RightTopOrientation:
      orient = mds_CreateScalar("righttop");
      break;

    case RightBottomOrientation:
      orient = mds_CreateScalar("rightbottom");
      break;

    case LeftBottomOrientation:
      orient = mds_CreateScalar("leftbottom");
      break;

    default:
      orient = mds_CreateScalar("undef");
      break;
  }
  E = ent_Create();
  ent_data (E) = orient;
  ent_type (E) = MATRIX_DENSE_STRING;
  install (bw, ("orientation"), E);

  // page
  size_t w0;
  size_t h0;
  ssize_t x0;
  ssize_t y0;
  MagickGetImagePage(magick_wand[default_wand_idx], &w0, &h0, &x0, &y0);
  MDR *pag  = mdi_Create(1,4);
  MdiV0(pag,0) = w0;
  MdiV0(pag,1) = h0;
  MdiV0(pag,2) = x0;
  MdiV0(pag,3) = y0;
  E = ent_Create();
  ent_data (E) = pag;
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("page"), E);

  // resolution
  MDR *resol  = mdr_Create(1,2);
  MagickGetImageResolution(magick_wand[default_wand_idx], &MdrV0(resol,0), &MdrV0(resol,1));
  E = ent_Create();
  ent_data (E) = resol;
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("resolution"), E);

  // for what follows:
  PixelWand *pw;
  MDR *rgba;

  // background
  pw = NewPixelWand();
  MagickGetImageBackgroundColor(magick_wand[default_wand_idx], pw);
  rgba = mdr_Create(1,4);
  MdrV0(rgba,0) = PixelGetRed     (pw);
  MdrV0(rgba,1) = PixelGetGreen   (pw);
  MdrV0(rgba,2) = PixelGetBlue    (pw);
  MdrV0(rgba,3) = PixelGetOpacity (pw);
  pw = DestroyPixelWand(pw);
  E = ent_Create();
  ent_data (E) = rgba;
  ent_type (E) = MATRIX_DENSE_REAL;
  install  (bw, "background", E);

  // border color
  pw = NewPixelWand();
  MagickGetImageBorderColor(magick_wand[default_wand_idx], pw);
  rgba = mdr_Create(1,4);
  MdrV0(rgba,0) = PixelGetRed     (pw);
  MdrV0(rgba,1) = PixelGetGreen   (pw);
  MdrV0(rgba,2) = PixelGetBlue    (pw);
  MdrV0(rgba,3) = PixelGetOpacity (pw);
  pw = DestroyPixelWand(pw);
  E = ent_Create();
  ent_data (E) = rgba;
  ent_type (E) = MATRIX_DENSE_REAL;
  install  (bw, "border_color", E);

  // image type
  ImageType imagetype=MagickGetImageType(magick_wand[default_wand_idx]);
  MDS *imtype=0;
  switch (imagetype)
  {
    case BilevelType:
      imtype = mds_CreateScalar("bilevel");
      break;

    case GrayscaleType:
      imtype = mds_CreateScalar("grayscale");
      break;

    case GrayscaleMatteType:
      imtype = mds_CreateScalar("grayscalematte");
      break;

    case PaletteType:
      imtype = mds_CreateScalar("palette");
      break;

    case PaletteMatteType:
      imtype = mds_CreateScalar("palettematte");
      break;

    case TrueColorType:
      imtype = mds_CreateScalar("truecolor");
      break;

    case TrueColorMatteType:
      imtype = mds_CreateScalar("truecolormatte");
      break;

    case ColorSeparationType:
      imtype = mds_CreateScalar("colorseparation");
      break;

    case ColorSeparationMatteType:
      imtype = mds_CreateScalar("colorseparationmatte");
      break;

    case OptimizeType:
      imtype = mds_CreateScalar("optimize");
      break;

    default:
      imtype = mds_CreateScalar("undef");
      break;
  }
  E = ent_Create();
  ent_data (E) = imtype;
  ent_type (E) = MATRIX_DENSE_STRING;
  install  (bw, "type", E);

  // units of resolution
  ResolutionType units=MagickGetImageUnits(magick_wand[default_wand_idx]);
  MDS *us=0;
  switch (units)
  {
    case PixelsPerInchResolution:
      us = mds_CreateScalar("ppi");
      break;

    case PixelsPerCentimeterResolution:
      us = mds_CreateScalar("ppcm");
      break;

    default:
      us = mds_CreateScalar("undef");
      break;
  }
  E = ent_Create();
  ent_data (E) = us;
  ent_type (E) = MATRIX_DENSE_STRING;
  install  (bw, "resolution_unit", E);

  // get image white point
  MDR *whp = mdr_Create(1,2);
  MagickGetImageWhitePoint(magick_wand[default_wand_idx], &MdrV0(whp,0),&MdrV0(whp,0));
  E = ent_Create();
  ent_data (E) = whp;
  ent_type (E) = MATRIX_DENSE_REAL;
  install  (bw, "white_point", E);


  rent = ent_Create ();
  ent_data (rent) = bw;
  ent_type (rent) = BTREE;
  return rent;
}

#if 0
//
// set default wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".draw"
#define IM_MAX_NO_METHOD_PARAMS 12
Ent * ent_im_wand_draw (int nargs, Datum args[])
{
  Ent *e1=0, *en=0, *rent;
  int i,j,rval=0;
  char *method=0;
  int np=0, have_ch=0;  // for counting parameters to the methods
  char *lb,*rb, *c=0;
  unsigned char ar=0,ac=0,ag=0,am=0,ab=0,ay=0,ao=0,aB=0,aM=0;
  double params[IM_MAX_NO_METHOD_PARAMS] = {0};
  MDR *mdr_param=0;

  // parameters needed for calling the functions
  void (*func)() = 0;      // function name
  MagickBooleanType (*func_ch)() = 0;   // function name
  NoiseType noise_type;                 // special function parameters

  MagickBooleanType bval;

  if (nargs < 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": ImageMagick wrapper library for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Draw method to active image in default Magick Wand.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(method,param1,....),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'method' is the name of the method to be applied, while.\n");
    fprintf (stdout,
             THIS_SOLVER ":   'param1,..' are the method parameters.\n");
    print_recognized_methods();
    fprintf (stdout,
             THIS_SOLVER ": Check manual for method parameters.\n");
    rerror ( THIS_SOLVER ": requires 1 or more arguments");
  }

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if (nargs >= 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_STRING)
      method  = class_char_pointer(e1);

    if (isvalidstring(method)<1)
      method = 0;
  }

  if (!method)
    rerror(THIS_SOLVER  ": First argument 'method' is required!\n");

  // identify the method and its arguments
  if (strstr(method,"clear"))
  {
    have_ch = 0;
    func = &ClearDrawingWand;
    np = -1;
  }
  if (strstr(method,"clone"))
  {
    have_ch = 0;
    func = &CloneDrawingWand;
    np = 0;
  }
  else if (strstr(method,"destroy"))
  {
    have_ch = 0;
    func = &DestroyDrawingWand;
    np = 0;
  }
  else if (strstr(method,"annote"))
  {
    // parameters: radius, sigma
    have_ch = 0;
    func = &DrawAnnotation;
    np = 0;
  }
  else if (strstr(method,"adaptivesharpen"))
  {
    // parameters: radius, sigma
    have_ch = 1;
    func = &MagickAdaptiveSharpenImage;
    func_ch = &MagickAdaptiveSharpenImageChannel;
    np = 2;
  }
  else if (strstr(method,"adaptivethreshold"))
  {
    // parameters: width, height, offset
    have_ch = 0;
    func = &_MagickAdaptiveThresholdImage;
    np = 3;
  }
  else if (strstr(method,"blackthreshold"))
  {
    // parameters: r,g,b,o
    have_ch = 0;
    func = &MagickBlackThresholdImage;
    np = -1;
  }
  else if (strstr(method,"threshold"))
  {
    // parameters: threshold
    have_ch = 1;
    func = &MagickThresholdImage;
    func_ch = &MagickThresholdImageChannel;
    np = 1;
  }
  else if (strstr(method,"addnoise"))
  {
    // parameters: "type", attanuate
    have_ch = 1;
    func = &MagickAddNoiseImage;
    func_ch = &MagickAddNoiseImageChannel;
    np = 2;
  }
  else if (strstr(method,"affine"))
  {
    // parameters: rx,sx,sy,ry,tx,ty
    have_ch = 0;
    func = &MagickAffineTransformImage;
    np = 6;
  }
  else if (strstr(method,"annotate"))
  {
    // parameters: text, [x,y], angle
    have_ch = 0;
    func = &MagickAnnotateImage;
    np = -1;
  }
  else if (strstr(method,"autogamma"))
  {
    have_ch = 1;
    func = &MagickAutoGammaImage;
    func_ch = &MagickAutoGammaImageChannel;
    np = 0;
  }
  else if (strstr(method,"autolevel"))
  {
    have_ch = 1;
    func = &MagickAutoLevelImage;
    func_ch = &MagickAutoLevelImageChannel;
    np = 0;
  }
  else if (strstr(method,"blueshift"))
  {
    // parameters: factor
    have_ch = 0;
    func = &MagickBlueShiftImage;
    np = 1;
  }
  else if (strstr(method,"gaussianblur"))
  {
    // parameters: radius,sigma
    have_ch = 1;
    func = &MagickGaussianBlurImage;
    func_ch = &MagickGaussianBlurImageChannel;
    np = 2;
  }
  else if (strstr(method,"motionblur"))
  {
    // parameters: radius, sigma, angle
    have_ch = 1;
    func = &MagickMotionBlurImage;
    func_ch = &MagickMotionBlurImageChannel;
    np = 3;
  }
  else if (strstr(method,"radialblur"))
  {
    // parameters: angle
    have_ch = 1;
    func = &MagickRadialBlurImage;
    func_ch = &MagickRadialBlurImageChannel;
    np = 1;
  }
  else if (strstr(method,"selectiveblur"))
  {
    // parameters: radius,sigma,threshold
    have_ch = 1;
    func = &MagickSelectiveBlurImage;
    func_ch = &MagickSelectiveBlurImageChannel;
    np = 3;
  }
  else if (strstr(method,"blur"))
  {
    // parameters: radius, sigma
    have_ch = 1;
    func = &MagickBlurImage;
    func_ch = &MagickBlurImageChannel;
    np = 2;
  }
  else if (strstr(method,"border"))
  {
    // parameters: r,g,b,o,w,h
    have_ch = 0;
    func = &MagickBorderImage;
    np = 6;
  }
  else if (strstr(method,"brightnesscontrast"))
  {
    // parameters: brightness, contrast
    have_ch = 1;
    func = &MagickBrightnessContrastImage;
    func_ch = &MagickBrightnessContrastImageChannel;
    np = 2;
  }
  else if (strstr(method,"charcoal"))
  {
    // parameters: radius, sigma
    have_ch = 0;
    func = &MagickCharcoalImage;
    np = 2;
  }
  else if (strstr(method,"chop"))
  {
    // parameters: width, height, x-offset, y-offset
    have_ch = 0;
    func = &MagickChopImage;
    np = 4;
  }
  else if (strstr(method,"cdl"))
  {
    // parameters: filename, or
    // offset[1x3], slope[1x3], power[1x3], saturation
    have_ch = 0;
    func = &MagickColorDecisionListImage;
    np = -1;
  }
  else if (strstr(method,"colorize"))
  {
    // parameters: filename, or
    // offset[1x3], slope[1x3], power[1x3], saturation
    have_ch = 0;
    func = &MagickColorizeImage;
    np = -1;
  }
  else if (strstr(method,"colormatrix"))
  {
    // parameters: filename, or
    // offset[1x3], slope[1x3], power[1x3], saturation
    have_ch = 0;
    func = &MagickColorMatrixImage;
    np = -1;
  }
  else if (strstr(method,"clamp"))
  {
    have_ch = 1;
    func = &MagickClampImage;
    func_ch = &MagickClampImageChannel;
    np = 0;
  }
  else if (strstr(method,"clip"))
  {
    have_ch = 0;
    func    = &MagickClipImage;
    np = 0;
  }
  else if (strstr(method,"contraststretch"))
  {
    // parameters: black,white
    have_ch = 1;
    func    = &MagickContrastStretchImage;
    func_ch = &MagickContrastStretchImageChannel;
    np = 2;
  }
  else if (strstr(method,"sigmoidalcontrast"))
  {
    // parameters: sharpen, alpha, beta
    have_ch = 1;
    func    = &_MagickSigmoidalContrastImage;
    func_ch = &_MagickSigmoidalContrastImageChannel;
    np = 3;
  }
  else if (strstr(method,"convolve"))
  {
    have_ch = 1;
    func    = &MagickConvolveImage;
    func    = &MagickConvolveImageChannel;
    np = -1;
  }
  else if (strstr(method,"contrast"))
  {
    have_ch = 0;
    func    = &_MagickContrastImage;
    np = 1;
  }
  else if (strstr(method,"crop"))
  {
    // parameters: width, height, x-offset, y-offset
    have_ch = 0;
    func = &MagickCropImage;
    np = 4;
  }
  else if (strstr(method,"cyclecolormap"))
  {
    // parameters: d
    have_ch = 0;
    func = &_MagickCycleColormapImage;
    np = 1;
  }
  else if (strstr(method,"deskew"))
  {
    // parameters: threshold
    have_ch = 0;
    func = &MagickDeskewImage;
    np = 1;
  }
  else if (strstr(method,"despecle"))
  {
    // parameters: threshold
    have_ch = 0;
    func = &MagickDespeckleImage;
    np = 0;
  }
  else if (strstr(method,"edge"))
  {
    // parameters: threshold
    have_ch = 0;
    func = &MagickEdgeImage;
    np = 1;
  }
  else if (strstr(method,"emboss"))
  {
    // parameters: radius, sigma
    have_ch = 0;
    func = &MagickEmbossImage;
    np = 2;
  }
  else if (strstr(method,"enhance"))
  {
    have_ch = 0;
    func = &MagickEnhanceImage;
    np = 0;
  }
  else if (strstr(method,"equalize"))
  {
    have_ch = 1;
    func = &MagickEqualizeImage;
    func_ch = &MagickEqualizeImageChannel;
    np = 0;
  }
  else if (strstr(method,"extent"))
  {
    // parameters: w,h,x,y
    have_ch = 0;
    func = &MagickExtentImage;
    np = 4;
  }
  else if (strstr(method,"flip"))
  {
    have_ch = 0;
    func = &MagickFlipImage;
    np = 0;
  }
  else if (strstr(method,"flop"))
  {
    have_ch = 0;
    func = &MagickFlopImage;
    np = 0;
  }
  else if (strstr(method,"gamma"))
  {
    // parameters: gamma
    have_ch = 1;
    func = &MagickGammaImage;
    func_ch = &MagickGammaImageChannel;
    np = 1;
  }
  else if (strstr(method,"implode"))
  {
    // parameters: gamma
    have_ch = 0;
    func = &MagickImplodeImage;
    np = 1;
  }
  else if (strstr(method,"level"))
  {
    // parameters: gamma
    have_ch = 1;
    func = &MagickLevelImage;
    func_ch = &MagickLevelImageChannel;
    np = 3;
  }
  else if (strstr(method,"linearstretch"))
  {
    // parameters: gamma
    have_ch = 0;
    func = &MagickLinearStretchImage;
    np = 2;
  }
  else if (strstr(method,"liquidrescale"))
  {
    // parameters: cols, rows, delta_x, rigidity
    have_ch = 0;
    func = &MagickLiquidRescaleImage;
    np = 4;
  }
  else if (strstr(method,"modulate"))
  {
    // parameters: brightness, saturation, hue
    have_ch = 0;
    func = &MagickModulateImage;
    np = 3;
  }
  else if (strstr(method,"negate"))
  {
    // parameters: gray
    have_ch = 1;
    func = &_MagickNegateImage;
    func_ch = &_MagickNegateImageChannel;
    np = 1;
  }
  else if (strstr(method,"normalize"))
  {
    have_ch = 1;
    func = &MagickNormalizeImage;
    func_ch = &MagickNormalizeImageChannel;
    np = 0;
  }
  else if (strstr(method,"oilpaint"))
  {
    // parameters: radius
    have_ch = 0;
    func = &MagickOilPaintImage;
    np = 1;
  }
  else if (strstr(method,"polaroid"))
  {
    // parameters: radius
    have_ch = 0;
    func = &_MagickPolaroidImage;
    np = 1;
  }
  else if (strstr(method,"posterize"))
  {
    // parameters: levels, dither
    have_ch = 0;
    func = &_MagickPosterizeImage;
    np = 2;
  }
  else if (strstr(method,"randomthreshold"))
  {
    // parameters: low, high
    have_ch = 1;
    func = &MagickRandomThresholdImage;
    func_ch = &MagickRandomThresholdImageChannel;
    np = 2;
  }
  else if (strstr(method,"resize"))
  {
    // parameters: low, high
    have_ch = 0;
    func = &MagickResizeImage;
    np = 4;
  }
  else if (strstr(method,"roll"))
  {
    // parameters: low, high
    have_ch = 0;
    func = &_MagickRollImage;
    np = 2;
  }
  else if (strstr(method,"rotate"))
  {
    // parameters: low, high
    have_ch = 0;
    func = &MagickRotateImage;
    np = -1;
  }
  else if (strstr(method,"sample"))
  {
    // parameters: low, high
    have_ch = 0;
    func = &_MagickSampleImage;
    np = 2;
  }
  else if (strstr(method,"scale"))
  {
    // parameters: low, high
    have_ch = 0;
    func = &_MagickScaleImage;
    np = 2;
  }
  else if (strstr(method,"sepia"))
  {
    // parameters: threshold
    have_ch = 0;
    func = &MagickSepiaToneImage;
    np = 1;
  }
  else if (strstr(method,"shade"))
  {
    // parameters: gray, azimuth, elevation
    have_ch = 0;
    func = &_MagickShadeImage;
    np = 3;
  }
  else if (strstr(method,"shadow"))
  {
    // parameters: opacity, sigma, x, y
    have_ch = 0;
    func = &_MagickShadowImage;
    np = 4;
  }
  else if (strstr(method,"shave"))
  {
    // parameters: opacity, sigma, x, y
    have_ch = 0;
    func = &_MagickShaveImage;
    np = 2;
  }
  else if (strstr(method,"shear"))
  {
    // parameters: opacity, sigma, x, y
    have_ch = 0;
    func = &MagickShearImage;
    np = -1;
  }
  else if (strstr(method,"sketch"))
  {
    // parameters: radius, sigma, angle
    have_ch = 0;
    func = &MagickSketchImage;
    np = 3;
  }
  else if (strstr(method,"solarize"))
  {
    // parameters: threshold
    have_ch = 0;
    func    = &MagickSolarizeImage;
    np = 1;
  }
  else if (strstr(method,"splice"))
  {
    // parameters: width, height, x, y
    have_ch = 0;
    func    = &_MagickSpliceImage;
    np = 4;
  }
  else if (strstr(method,"spread"))
  {
    // parameters: radius
    have_ch = 0;
    func    = &MagickSpreadImage;
    np = 1;
  }
  else if (strstr(method,"strip"))
  {
    have_ch = 0;
    func    = &MagickStripImage;
    np = 0;
  }
  else if (strstr(method,"swirl"))
  {
    // parameters: radius
    have_ch = 0;
    func    = &MagickSwirlImage;
    np = 1;
  }
  else if (strstr(method,"thumbnail"))
  {
    // parameters: radius
    have_ch = 0;
    func    = &_MagickThumbnailImage;
    np = 2;
  }
  else if (strstr(method,"transparentpaint"))
  {
    have_ch = 0;
    func    = &MagickTransparentPaintImage;
    np = -1;
  }
  else if (strstr(method,"transpose"))
  {
    have_ch = 0;
    func    = &MagickTransposeImage;
    np = 0;
  }
  else if (strstr(method,"transverse"))
  {
    have_ch = 0;
    func    = &MagickTransverseImage;
    np = 0;
  }
  else if (strstr(method,"trim"))
  {
    // parameters: fuzz
    have_ch = 0;
    func    = &MagickTrimImage;
    np = 1;
  }
  else if (strstr(method,"uniquecolors"))
  {
    have_ch = 0;
    func    = &MagickUniqueImageColors;
    np = 0;
  }
  else if (strstr(method,"unsharpmask"))
  {
    // parameters: radius, sigma, amount, threshold
    have_ch = 1;
    func    = &MagickUnsharpMaskImage;
    func_ch = &MagickUnsharpMaskImageChannel;
    np = 4;
  }
  else if (strstr(method,"vignette"))
  {
    have_ch = 0;
    func = &_MagickVignetteImage;
    np = 4;
  }
  else if (strstr(method,"wave"))
  {
    have_ch = 0;
    func = &MagickWaveImage;
    np = 2;
  }
  else if (strstr(method,"whitethreshold"))
  {
    have_ch = 0;
    func = &MagickWhiteThresholdImage;
    np = -1;
  }
  else
  {
    fprintf(stderr, "Method %s not recognized!\n", method);
    print_recognized_methods();
    rerror(THIS_SOLVER  ": Cannot continue. Method is missing!\n");
  }

  //
  // process channels if necessary
  //
  if (have_ch)
  {
    lb = strstr(method,"[");
    rb = strstr(method,"]");
    if (lb && rb)
    {
      ar = (strstr(lb, "r") != 0);  // RedChannel
      ac = (strstr(lb, "c") != 0);  // CyanChannel
      ag = (strstr(lb, "g") != 0);  // GreenChannel
      am = (strstr(lb, "m") != 0);  // MagentaChannel
      ab = (strstr(lb, "b") != 0);  // BlueChannel
      ay = (strstr(lb, "y") != 0);  // YellowChannel
      ao = (strstr(lb, "o") != 0);  // OpacityChannel
      aB = (strstr(lb, "B") != 0);  // BlackChannel
      aM = (strstr(lb, "M") != 0);  // MatteChannel
    }
  }

  //
  // first process functions with special arguments
  // (that cannot be casted to double * for arguments)
  //
  if(func == &MagickAddNoiseImage)
  {
    if(nargs != np+1 && nargs != np)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_STRING)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires string as a second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    c = class_char_pointer(en);
    switch(*c)
    {
      case 'u':
      case 'U':
        noise_type = UniformNoise;
        break;

      case 'g':
      case 'G':
      case 'n':
      case 'N':
        noise_type = GaussianNoise;
        break;

      case 'm':
      case 'M':
        noise_type = MultiplicativeGaussianNoise;
        break;

      case 'i':
      case 'I':
        noise_type = ImpulseNoise;
        break;

      case 'l':
      case 'L':
        noise_type = LaplacianNoise;
        break;

      case 'p':
      case 'P':
        noise_type = PoissonNoise;
        break;

      default:
        noise_type = UniformNoise;
        break;
    }
    ent_Clean(en);

    if (nargs==3)
    {
      en = bltin_get_ent (args[2]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      sprintf(default_text, "%g", class_double(en));
      ent_Clean(en);

      MagickSetImageArtifact(magick_wand[default_wand_idx], "attenuate", default_text);
    }

    if (ar)
      func_ch(magick_wand[default_wand_idx], RedChannel, noise_type);
    if (ac)
      func_ch(magick_wand[default_wand_idx], CyanChannel,noise_type);
    if (ag)
      func_ch(magick_wand[default_wand_idx], GreenChannel,noise_type);
    if (am)
      func_ch(magick_wand[default_wand_idx], MagentaChannel,noise_type);
    if (ab)
      func_ch(magick_wand[default_wand_idx], BlueChannel,noise_type);
    if (ay)
      func_ch(magick_wand[default_wand_idx], YellowChannel,noise_type);
    if (ao)
      func_ch(magick_wand[default_wand_idx], OpacityChannel,noise_type);
    if (aB)
      func_ch(magick_wand[default_wand_idx], BlackChannel,noise_type);
    if (aM)
      func_ch(magick_wand[default_wand_idx], MatteChannel,noise_type);
    if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
      func(magick_wand[default_wand_idx], noise_type);
  }
  else if (func == &MagickResizeImage)
  {
    // special methods:
    if(nargs != 4)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x2 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires two entries vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=2)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires two entries vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = MdrV0(mdr_param,0);
    params[1] = MdrV0(mdr_param,1);
    ent_Clean(en);

    // third parameter is the filter name
    en = bltin_get_ent (args[2]);
    if (ent_type(en) != MATRIX_DENSE_STRING)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires string as a second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    FilterTypes filter;
    c = class_char_pointer(en);
    if (strstr(c, "bes"))
      filter = BesselFilter;
    else if (strstr(c, "bl"))
      filter = BlackmanFilter;
    else if (strstr(c,"box"))
      filter = BoxFilter;
    else if (strstr(c,"cat"))
      filter = CatromFilter;
    else if (strstr(c,"cub"))
      filter = CubicFilter;
    else if (strstr(c,"han"))
      filter = HanningFilter;
    else if (strstr(c,"her"))
      filter = HermiteFilter;
    else if (strstr(c,"lan"))
      filter = LanczosFilter;
    else if (strstr(c,"mit"))
      filter = MitchellFilter;
    else if (strstr(c,"quad"))
      filter = QuadraticFilter;
    else if (strstr(c,"sinc"))
      filter = SincFilter;
    else if (strstr(c,"tri"))
      filter = TriangleFilter;
    else
      filter = LanczosFilter;

    ent_Clean(en);

    // fourth is the blur
    en = bltin_get_ent (args[3]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires real 'blur' as fourth argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[2] = class_double(en);
    ent_Clean(en);

    func(magick_wand[default_wand_idx], (size_t) params[0], (size_t) params[1], filter, params[2]);
  }
  else if (func == &MagickRotateImage)
  {
    if(nargs != 3)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=3 && mdr_param->nrow * mdr_param->ncol !=4)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = MdrV0(mdr_param,0);
    params[1] = MdrV0(mdr_param,1);
    params[2] = MdrV0(mdr_param,2);
    if (mdr_param->nrow * mdr_param->ncol ==4 )
      params[3] = MdrV0(mdr_param,3);
    else
      params[3] = 1;
    ent_Clean(en);

    // fourth is the blur
    en = bltin_get_ent (args[2]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[4] = class_double(en);
    ent_Clean(en);

    PixelWand *pw = NewPixelWand();
    PixelSetRed(pw, params[0]);
    PixelSetGreen(pw, params[1]);
    PixelSetBlue(pw, params[2]);
    PixelSetOpacity(pw, params[3]);
    MagickRotateImage(magick_wand[default_wand_idx], pw,params[4]);
    pw = DestroyPixelWand(pw);
  }
  else if (func == &MagickShearImage)
  {
    if(nargs != 4)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=3 && mdr_param->nrow * mdr_param->ncol !=4)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = MdrV0(mdr_param,0);
    params[1] = MdrV0(mdr_param,1);
    params[2] = MdrV0(mdr_param,2);
    if (mdr_param->nrow * mdr_param->ncol ==4 )
      params[3] = MdrV0(mdr_param,3);
    else
      params[3] = 1;
    ent_Clean(en);

    // third is x_shear
    en = bltin_get_ent (args[2]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[4] = class_double(en);
    ent_Clean(en);

    // fourth is y_shear
    en = bltin_get_ent (args[3]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[5] = class_double(en);
    ent_Clean(en);

    PixelWand *pw = NewPixelWand();
    PixelSetRed(pw, params[0]);
    PixelSetGreen(pw, params[1]);
    PixelSetBlue(pw, params[2]);
    PixelSetOpacity(pw, params[3]);

    MagickShearImage(magick_wand[default_wand_idx], pw, params[4], params[5]);

    pw = DestroyPixelWand(pw);
  }
  else if (func == &MagickTransparentPaintImage)
  {
    if(nargs != 5)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
      en = bltin_get_ent (args[1]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      mdr_param = ent_data(en);
      if (mdr_param->nrow * mdr_param->ncol !=3 && mdr_param->nrow * mdr_param->ncol !=4)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[0] = MdrV0(mdr_param,0);
      params[1] = MdrV0(mdr_param,1);
      params[2] = MdrV0(mdr_param,2);
      if (mdr_param->nrow * mdr_param->ncol ==4 )
        params[3] = MdrV0(mdr_param,3);
      else
        params[3] = 1;
      ent_Clean(en);

      // third is alpha
      en = bltin_get_ent (args[2]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[4] = class_double(en);
      ent_Clean(en);

      // fourth is fuzz
      en = bltin_get_ent (args[3]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[5] = class_double(en);
      ent_Clean(en);

      // fifth is invert
      en = bltin_get_ent (args[4]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires real 'angle' as third argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[6] = class_double(en);
      ent_Clean(en);

      PixelWand *pw = NewPixelWand();
      PixelSetRed(pw, params[0]);
      PixelSetGreen(pw, params[1]);
      PixelSetBlue(pw, params[2]);
      PixelSetOpacity(pw, params[3]);

      MagickTransparentPaintImage(magick_wand[default_wand_idx], pw, params[4], params[5],
                                  (MagickBooleanType) params[6]);

      pw = DestroyPixelWand(pw);
  }
  else if (func == &MagickWhiteThresholdImage)
  {
    if(nargs != 2)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=3 && mdr_param->nrow * mdr_param->ncol !=4)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = MdrV0(mdr_param,0);
    params[1] = MdrV0(mdr_param,1);
    params[2] = MdrV0(mdr_param,2);
    if (mdr_param->nrow * mdr_param->ncol ==4 )
      params[3] = MdrV0(mdr_param,3);
    else
      params[3] = 1;
    ent_Clean(en);

    PixelWand *pw = NewPixelWand();
    PixelSetRed(pw, params[0]);
    PixelSetGreen(pw, params[1]);
    PixelSetBlue(pw, params[2]);
    PixelSetOpacity(pw, params[3]);

    MagickWhiteThresholdImage(magick_wand[default_wand_idx], pw);

    pw = DestroyPixelWand(pw);
  }
  else if (func == &MagickBlackThresholdImage)
  {
    if(nargs != 2)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=3 && mdr_param->nrow * mdr_param->ncol !=4)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] real vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = MdrV0(mdr_param,0);
    params[1] = MdrV0(mdr_param,1);
    params[2] = MdrV0(mdr_param,2);
    if (mdr_param->nrow * mdr_param->ncol ==4 )
      params[3] = MdrV0(mdr_param,3);
    else
      params[3] = 1;
    ent_Clean(en);

    PixelWand *pw = NewPixelWand();
    PixelSetRedQuantum(pw, (Quantum) params[0]);
    PixelSetGreenQuantum(pw, (Quantum) params[1]);
    PixelSetBlueQuantum(pw, (Quantum) params[2]);
    PixelSetOpacityQuantum(pw, (Quantum) params[3]);

    MagickBlackThresholdImage(magick_wand[default_wand_idx], pw);

    pw = DestroyPixelWand(pw);
  }
  else if (func == &MagickAnnotateImage)
  {
    if(nargs != 4)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // first argument: text
    char *ann=0;
    Ent *es = bltin_get_ent (args[1]);
    if (ent_type(es) != MATRIX_DENSE_STRING)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires text as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    ann = class_char_pointer(es);

    // second argument: [x,y]
    en = bltin_get_ent (args[2]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires [x,y] real vector as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow * mdr_param->ncol !=2)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires [x,y] real vector as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = MdrV0(mdr_param,0);
    params[1] = MdrV0(mdr_param,1);
    ent_Clean(en);

    // third argument: angle
    en = bltin_get_ent (args[3]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires real scalar 'angle' as fourth argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[2] = class_double(en);
    ent_Clean(en);

    DrawingWand *drawing_wand = NewDrawingWand();
    rval = MagickAnnotateImage(magick_wand[default_wand_idx],drawing_wand,
                               params[0],params[1],params[2],ann);
    drawing_wand = DestroyDrawingWand(drawing_wand);
    ent_Clean(es);
  }
  else if (func == &MagickColorDecisionListImage)
  {
    char *filename=0;

    if(nargs != 2 && nargs!=5)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    if (nargs == 2)
    {
      Ent *en = bltin_get_ent (args[1]);
      if (ent_type(en) != MATRIX_DENSE_STRING)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires filename as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      filename = class_char_pointer(en);
      if (!filename)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires filename as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      FILE *fr = fopen (filename, "r");
      if (!fr)
      {
        fprintf(stderr, THIS_SOLVER  ": file %s does not exist!\n", filename);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      i=0;
      while( ( default_text[i] = fgetc(fr) ) != EOF  && i<MAX_LEN_DEFAULT_TEXT)
        i++;
      if (i==MAX_LEN_DEFAULT_TEXT)
      {
        fprintf(stderr, THIS_SOLVER  ": file %s is too big! It should contain a single profile!\n", filename);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      default_text[i]='\0';
      fclose(fr);
      ent_Clean(en);
    }
    else
    {
      // we create XML string from the user provided input data

      // second argument is a 1x3 matrix
      en = bltin_get_ent (args[1]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'offset' real vector as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      mdr_param = ent_data(en);
      if (mdr_param->nrow * mdr_param->ncol !=3)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'offset' real vector as second argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[0] = MdrV0(mdr_param,0);
      params[1] = MdrV0(mdr_param,1);
      params[2] = MdrV0(mdr_param,2);
      ent_Clean(en);

      // third argument is a 1x3 matrix
      en = bltin_get_ent (args[2]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'slope' real vector as third argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      mdr_param = ent_data(en);
      if (mdr_param->nrow * mdr_param->ncol !=3)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'slope' real vector as third argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[3] = MdrV0(mdr_param,0);
      params[4] = MdrV0(mdr_param,1);
      params[5] = MdrV0(mdr_param,2);
      ent_Clean(en);

      // fourth argument is a 1x3 matrix
      en = bltin_get_ent (args[3]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'power' real vector as fourth argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      mdr_param = ent_data(en);
      if (mdr_param->nrow * mdr_param->ncol !=3)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'power' real vector as fourth argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[6] = MdrV0(mdr_param,0);
      params[7] = MdrV0(mdr_param,1);
      params[8] = MdrV0(mdr_param,2);
      ent_Clean(en);

      // fifth argument is a scalar
      en = bltin_get_ent (args[4]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
      {
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'saturation' real scalar as fifth argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[9] = class_double(en);
      ent_Clean(en);

      // now assemble XML file
      i=0;
      char *current_text=default_text;
      //
      sprintf(&current_text[i],"<ColorCorrectionCollection xmlns=\"urn:ASC:CDL:v0.01\">\n");
      i = strlen(current_text);
      //
      sprintf(&current_text[i],"<ColorCorrection id=\"rlab0000\">\n");
      i = strlen(current_text);
      //
      sprintf(&current_text[i],"<SOPNode>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"<Description>RLaB was here!</Description>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"<Offset> %g %g %g </Offset>\n",params[0],params[1],params[2]);
      i = strlen(current_text);
      sprintf(&current_text[i],"<Slope> %g %g %g </Slope>\n",params[3],params[4],params[5]);
      i = strlen(current_text);
      sprintf(&current_text[i],"<Power> %g %g %g </Power>\n",params[6],params[7],params[8]);
      i = strlen(current_text);
      sprintf(&current_text[i],"</SOPNode>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"<SatNode>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"<Saturation> %g </Saturation>\n",params[9]);
      i = strlen(current_text);
      sprintf(&current_text[i],"</SatNode>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"</ColorCorrection>\n");
      i = strlen(current_text);
      sprintf(&current_text[i],"</ColorCorrectionCollection>\n");
      i = strlen(current_text);
      current_text[i] = '\0';
    }

    if (im_debug)
      fprintf(stderr, "%s", default_text);

    bval = MagickColorDecisionListImage(magick_wand[default_wand_idx], default_text);
    if (bval == MagickFalse)
      WriteWandException(magick_wand[default_wand_idx]);
  }
  else if (func == &MagickColorizeImage)
  {
    if(nargs != 2)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 or 1x4 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'color' real 3- or 4-vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow *mdr_param->ncol != 3 && mdr_param->nrow *mdr_param->ncol != 4)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'color' real 3- or 4-vector as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }

    PixelWand *c = NewPixelWand();
    PixelSetRed   (c, MdrV0(mdr_param,0));
    PixelSetGreen (c, MdrV0(mdr_param,1));
    PixelSetBlue  (c, MdrV0(mdr_param,2));

    PixelWand *o = NewPixelWand();
    if (mdr_param->nrow *mdr_param->ncol == 4)
      PixelSetOpacity(o, MdrV0(mdr_param,3));
    else
      PixelSetOpacity(o, 1.0);

    MagickColorizeImage(magick_wand[default_wand_idx], c, o);

    c = DestroyPixelWand(c);
    o = DestroyPixelWand(o);
  }
  else if (func == &MagickColorMatrixImage)
  {
    if(nargs != 2)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    // second argument is a 1x3 matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow != mdr_param->ncol)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    if (mdr_param->nrow != 3 && mdr_param->nrow != 4 && mdr_param->nrow != 5)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }

    // construct string: rows in the matrix
    i=0;
    char *current_text=default_text;

    for (i=0; i<mdr_param->nrow; i++)
      for (j=0; j<mdr_param->ncol; j++)
      {
        sprintf(&current_text[i],"%g ", Mdr0(mdr_param,i,j));
        i = strlen(current_text);
      }
    current_text[i-1] ='\0';  // remove extraspace at the end

    KernelInfo *mykernel = AcquireKernelInfo(default_text);
    MagickColorMatrixImage(magick_wand[default_wand_idx],mykernel);

    mykernel = DestroyKernelInfo(mykernel);
    ent_Clean(en);
  }
  else if  (func == &MagickConvolveImage)
  {
    if(nargs != 2)
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    MDR *order=0;

    // second argument is a square matrix
    en = bltin_get_ent (args[1]);
    if (ent_type(en) != MATRIX_DENSE_REAL)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow != mdr_param->ncol)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    order = mdr_Transpose(mdr_param);

    if (ar)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], RedChannel, (size_t) MNR(order), MDRPTR(order));
    if (ac)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], CyanChannel, (size_t) MNR(order), MDRPTR(order));
    if (ag)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], GreenChannel, (size_t) MNR(order), MDRPTR(order));
    if (am)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], MagentaChannel, (size_t) MNR(order), MDRPTR(order));
    if (ab)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], BlueChannel, (size_t) MNR(order), MDRPTR(order));
    if (ay)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], YellowChannel, (size_t) MNR(order), MDRPTR(order));
    if (ao)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], OpacityChannel, (size_t) MNR(order), MDRPTR(order));
    if (aB)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], BlackChannel, (size_t) MNR(order), MDRPTR(order));
    if (aM)
      MagickConvolveImageChannel(magick_wand[default_wand_idx], MatteChannel, (size_t) MNR(order), MDRPTR(order));
    if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
      MagickConvolveImage(magick_wand[default_wand_idx], (size_t) MNR(order), MDRPTR(order));

    mdr_Destroy(order);
    ent_Clean(en);
  }
  else
  {
    //
    // generic method with generic parameters
    //
    if (nargs==2 && np>1)
    {
      // all parameters are in the single array
      en = bltin_get_ent (args[1]);
      if (ent_type(en) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

      mdr_param = ent_data(en);
      if (mdr_param->nrow * mdr_param->ncol!=np)
        rerror(THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      for (i=0; i<np; i++)
        params[i] = MdrV0(mdr_param, i);

      ent_Clean(en);
    }
    else if (nargs==(np+1) && np>1)
    {
      // read parameters from arguments to the function
      for (i=1; i<=np; i++)
      {
        en = bltin_get_ent (args[i]);
        if (ent_type(en) == MATRIX_DENSE_REAL)
          params[i-1] = class_double(en);
        ent_Clean(en);
      }
    }
    else if ((np+1)!=nargs)
      rerror(THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");

    //
    // generic method with generic parameters
    //
    switch (np)
    {
      case 0:
        // no parameter functions
        if (!have_ch)
        {
          func(magick_wand[default_wand_idx]);
        }
        else
        {
          if (ar)
            func_ch(magick_wand[default_wand_idx], RedChannel);
          if (ac)
            func_ch(magick_wand[default_wand_idx], CyanChannel);
          if (ag)
            func_ch(magick_wand[default_wand_idx], GreenChannel);
          if (am)
            func_ch(magick_wand[default_wand_idx], MagentaChannel);
          if (ab)
            func_ch(magick_wand[default_wand_idx], BlueChannel);
          if (ay)
            func_ch(magick_wand[default_wand_idx], YellowChannel);
          if (ao)
            func_ch(magick_wand[default_wand_idx], OpacityChannel);
          if (aB)
            func_ch(magick_wand[default_wand_idx], BlackChannel);
          if (aM)
            func_ch(magick_wand[default_wand_idx], MatteChannel);
          if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
            func(magick_wand[default_wand_idx]);
        }
        break;

      case 1:
        if (!have_ch)
        {
          func(magick_wand[default_wand_idx], params[0]);
        }
        else
        {
          if (ar)
            func_ch(magick_wand[default_wand_idx], RedChannel, params[0]);
          if (ac)
            func_ch(magick_wand[default_wand_idx], CyanChannel,params[0]);
          if (ag)
            func_ch(magick_wand[default_wand_idx], GreenChannel,params[0]);
          if (am)
            func_ch(magick_wand[default_wand_idx], MagentaChannel,params[0]);
          if (ab)
            func_ch(magick_wand[default_wand_idx], BlueChannel,params[0]);
          if (ay)
            func_ch(magick_wand[default_wand_idx], YellowChannel,params[0]);
          if (ao)
            func_ch(magick_wand[default_wand_idx], OpacityChannel,params[0]);
          if (aB)
            func_ch(magick_wand[default_wand_idx], BlackChannel,params[0]);
          if (aM)
            func_ch(magick_wand[default_wand_idx], MatteChannel,params[0]);
          if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
            func(magick_wand[default_wand_idx], params[0]);
        }
        break;

      case 2:
        // generic two parameter functions
        if (!have_ch)
        {
          func(magick_wand[default_wand_idx], params[0], params[1]);
        }
        else
        {
          if (ar)
            func_ch(magick_wand[default_wand_idx], RedChannel,params[0], params[1]);
          if (ac)
            func_ch(magick_wand[default_wand_idx], CyanChannel,params[0], params[1]);
          if (ag)
            func_ch(magick_wand[default_wand_idx], GreenChannel,params[0], params[1]);
          if (am)
            func_ch(magick_wand[default_wand_idx], MagentaChannel,params[0], params[1]);
          if (ab)
            func_ch(magick_wand[default_wand_idx], BlueChannel,params[0], params[1]);
          if (ay)
            func_ch(magick_wand[default_wand_idx], YellowChannel,params[0], params[1]);
          if (ao)
            func_ch(magick_wand[default_wand_idx], OpacityChannel,params[0], params[1]);
          if (aB)
            func_ch(magick_wand[default_wand_idx], BlackChannel,params[0], params[1]);
          if (aM)
            func_ch(magick_wand[default_wand_idx], MatteChannel,params[0], params[1]);
          if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
            func(magick_wand[default_wand_idx], params[0], params[1]);
        }
        break;

      case 3:
        // three parameter functions
        if (!have_ch)
          func(magick_wand[default_wand_idx], params[0], params[1],params[2]);
        else
        {
          if (ar)
            func_ch(magick_wand[default_wand_idx],RedChannel,params[0],params[1],params[2]);
          if (ac)
            func_ch(magick_wand[default_wand_idx],CyanChannel,params[0],params[1],params[2]);
          if (ag)
            func_ch(magick_wand[default_wand_idx],GreenChannel,params[0],params[1],params[2]);
          if (am)
            func_ch(magick_wand[default_wand_idx],MagentaChannel,params[0],params[1],params[2]);
          if (ab)
            func_ch(magick_wand[default_wand_idx],BlueChannel,params[0],params[1],params[2]);
          if (ay)
            func_ch(magick_wand[default_wand_idx],YellowChannel,params[0],params[1],params[2]);
          if (ao)
            func_ch(magick_wand[default_wand_idx],OpacityChannel,params[0],params[1],params[2]);
          if (aB)
            func_ch(magick_wand[default_wand_idx],BlackChannel,params[0],params[1],params[2]);
          if (aM)
            func_ch(magick_wand[default_wand_idx],MatteChannel,params[0],params[1],params[2]);
          if(!ar && !ac && !ag && !am && !ab && !ay && !ao && !aB && !aM)
            func(magick_wand[default_wand_idx], params[0], params[1],params[2]);
        }
        break;

      case 4:
        if (func == &MagickChopImage || func == &MagickCropImage || func == &MagickExtentImage)
        {
          func(magick_wand[default_wand_idx], (size_t) params[0], (size_t) params[1],
               (size_t) params[2], (size_t) params[3]);
        }
        else if (func == &MagickLiquidRescaleImage)
        {
          func(magick_wand[default_wand_idx], (size_t) params[0], (size_t) params[1], params[2],params[3]);
        }
        break;

      case 5:
        break;

      case 6:
        if (func == &MagickAffineTransformImage)
        {
          fprintf(stderr, "'affinetransform' is not working!\n");
          AffineMatrix aff_mat;
          DrawingWand *aff_dw=NewDrawingWand();
          aff_mat.sx = params[0];
          aff_mat.rx = params[1];
          aff_mat.ry = params[2];
          aff_mat.sy = params[3];
          aff_mat.tx = params[4];
          aff_mat.ty = params[5];
          DrawAffine(aff_dw, &aff_mat);
//           DrawRender(aff_dw);
          if (func(magick_wand[default_wand_idx], aff_dw) == MagickFalse)
            WriteWandException(magick_wand[default_wand_idx]);
//           DrawRender(aff_dw);
//           MagickDrawImage(magick_wand[default_wand_idx],aff_dw);
          aff_dw = DestroyDrawingWand(aff_dw);
        }
        else if (func == &MagickBorderImage)
        {
          PixelWand *pw = NewPixelWand();
          PixelSetRed     (pw, params[0]);
          PixelSetGreen   (pw, params[1]);
          PixelSetBlue    (pw, params[2]);
          PixelSetOpacity (pw, params[3]);
          func(magick_wand[default_wand_idx], pw, (size_t) params[4], (size_t) params[5]);
          pw = DestroyPixelWand(pw);
        }
        break;

      default:
        break;
    }
  }

  ent_Clean(e1);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (rval);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

#endif

#endif  /* HAVE_IMAGEMAGICK */
