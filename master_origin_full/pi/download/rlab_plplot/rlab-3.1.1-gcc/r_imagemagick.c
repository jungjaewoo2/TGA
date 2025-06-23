//
// rimagemagick.c:  rlab's interface to ImageMagick 6 and 7
// Marijan Kostrun, 2013-2018
//
// This file is a part of RLaB + rlabplus
// Copyright (C) 2018  Marijan Kostrun
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
#include "buffer.h"

#include <stdio.h>
#include <fcntl.h>

//
#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"
#include "rlab_macros_code.h"

#ifdef HAVE_IMAGEMAGICK  /* Set in config.h */

#include <wand/MagickWand.h>

#define MAX_NUMBER_MAGICK_WANDS 32
static int InitMagickWandLib=0;
static int default_wand_idx=-1;
static MagickWand *magick_wand[MAX_NUMBER_MAGICK_WANDS] = {0};
static unsigned char im_debug=1;
static ExceptionType severity;
static char *description=0;

#define RLAB_MAGICK_QUANTUM_SIZE 255

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


static const char * imagetype2string(ImageType imval)
{
  switch(imval)
  {
    case BilevelType:
      return "bilevel";

    case GrayscaleType:
      return "grayscale";

    case PaletteType:
      return "palette";

    case TrueColorType:
      return "truecolor";

    case ColorSeparationType:
      return "colorseparation";

    case OptimizeType:
      return "optimize";

#if (HAVE_IMAGEMAGICK == 6)

    case GrayscaleMatteType:
      return "grayscalematte";

    case PaletteMatteType:
      return "palettematte";

    case TrueColorMatteType:
      return "truecolormatte";

    case ColorSeparationMatteType:
      return "colorseparationmatte";

    case PaletteBilevelMatteType:
      return "palettebilevelmatte";

#else

    case GrayscaleAlphaType:
      return "grayscalealpha";

    case PaletteAlphaType:
      return "palettealpha";

    case TrueColorAlphaType:
      return "truecoloralpha";

    case ColorSeparationAlphaType:
      return "colorseparationalpha";

    case PaletteBilevelAlphaType:
      return "palettebilevelalpha";

#endif

    default:
      return "undef";
  }

  return "undef";
}


static ImageType string2imagetype(char *val)
{
  if (isvalidstring(val) < 1)
    return UndefinedType;

  if (!strcmp(val,"bilevel"))
    return BilevelType;

  if (!strcmp(val,"grayscale"))
    return GrayscaleType;

  if (!strcmp(val,"palette"))
    return PaletteType;

  if (!strcmp(val,"truecolor"))
    return TrueColorType;

  if (!strcmp(val,"colorseparation"))
    return ColorSeparationType;

  if (!strcmp(val,"optimize"))
    return OptimizeType;

#if (HAVE_IMAGEMAGICK == 6)

  if (!strcmp(val,"grayscalematte"))
    return GrayscaleMatteType;

  if (!strcmp(val,"palettematte"))
    return PaletteMatteType;

  if (!strcmp(val,"truecolormatte"))
    return TrueColorMatteType;

  if (!strcmp(val,"colorseparationmatte"))
    return ColorSeparationMatteType;

  if (!strcmp(val,"palettebilevelmatte"))
    return PaletteBilevelMatteType;

#else

  if (!strcmp(val,"grayscalealpha"))
    return GrayscaleAlphaType;

  if (!strcmp(val,"palettealpha"))
    return PaletteAlphaType;

  if (!strcmp(val,"truecoloralpha"))
    return TrueColorAlphaType;

  if (!strcmp(val,"colorseparationalpha"))
    return ColorSeparationAlphaType;

  if (!strcmp(val,"palettebilevelalpha"))
    return PaletteBilevelAlphaType;

#endif

  return UndefinedType;
}


static ColorspaceType string2colorspace(char * val)
{
  if (isvalidstring(val)<1)
    return UndefinedColorspace;

  if (!strcmp(val,"cmyk"))
    return CMYKColorspace;

  if (!strcmp(val,"cmy"))
    return CMYColorspace;

  if (!strcmp(val,"gray"))
    return GRAYColorspace;

   if (!strcmp(val,"hclp"))
     return HCLpColorspace;

   if (!strcmp(val,"hcl"))
     return HCLColorspace;

   if (!strcmp(val,"hsb"))
     return HSBColorspace;

   if (!strcmp(val,"hsi"))
     return HSIColorspace;

   if (!strcmp(val,"hsl"))
     return HSLColorspace;

   if (!strcmp(val,"hsv"))
     return HSVColorspace;

   if (!strcmp(val,"hwb"))
     return HWBColorspace;

   if (!strcmp(val,"lab"))
     return LabColorspace;

   if (!strcmp(val,"lchab"))
     return LCHabColorspace;

   if (!strcmp(val,"lchuv"))
     return LCHuvColorspace;

   if (!strcmp(val,"lch"))
     return LCHColorspace;

   if (!strcmp(val,"log"))
     return LogColorspace;

   if (!strcmp(val,"lms"))
     return LMSColorspace;

   if (!strcmp(val,"luv"))
     return LuvColorspace;

   if (!strcmp(val,"ohta"))
    return OHTAColorspace;

   if (!strcmp(val,"xyz"))
    return XYZColorspace;

   if (!strcmp(val,"rec601ycbcr"))
     return Rec601YCbCrColorspace;

   if (!strcmp(val,"rec709ycbcr"))
     return Rec709YCbCrColorspace;

#if (HAVE_IMAGEMAGICK == 6)

   if (!strcmp(val,"rec601luma"))
     return Rec601LumaColorspace;

   if (!strcmp(val,"rec709luma"))
     return Rec709LumaColorspace;

#else

   if (!strcmp(val,"lgray"))
     return LinearGRAYColorspace;

   if (!strcmp(val,"xyy"))
     return xyYColorspace;

#endif

   if (!strcmp(val,"scrgb"))
     return scRGBColorspace;

   if (!strcmp(val,"srgb"))
     return sRGBColorspace;

   if (!strcmp(val,"rgb" ))
     return RGBColorspace;

   if (!strcmp(val,"transparent"))
     return TransparentColorspace;

   if (!strcmp(val,"ycbcr"))
    return YCbCrColorspace;

   if (!strcmp(val,"ycc"))
    return YCCColorspace;

   if (!strcmp(val,"ydbdr"))
     return YDbDrColorspace;

   if (!strcmp(val,"yiq"))
    return YIQColorspace;

   if (!strcmp(val,"ypbpr"))
    return YPbPrColorspace;

   if (!strcmp(val,"yuv"))
     return YUVColorspace;

   return UndefinedColorspace;
}

//
//
//
static const char * colorspace2string(ColorspaceType cs)
{
  switch (cs)
  {
    case CMYKColorspace:
      return "cmyk";

    case CMYColorspace:
      return "cmy";

    case GRAYColorspace:
      return "gray";

    case HCLpColorspace:
      return "hclp";

    case HCLColorspace:
      return "hcl";

    case HSBColorspace:
      return "hsb";

    case HSIColorspace:
      return "hsi";

    case HSLColorspace:
      return "hsl";

    case HSVColorspace:
      return "hsv";

    case HWBColorspace:
      return "hwb";

    case LabColorspace:
      return "lab";

    case LCHabColorspace:
      return "lchab";

    case LCHuvColorspace:
      return "lchuv";

    case LCHColorspace:
      return "lch";

    case LogColorspace:
      return "log";

    case LMSColorspace:
      return "lms";

    case LuvColorspace:
      return "luv";

    case OHTAColorspace:
      return "ohta";

    case XYZColorspace:
      return "xyz";

    case Rec601YCbCrColorspace:
      return "rec601ycbcr";

    case Rec709YCbCrColorspace:
      return "rec709ycbcr";

#if (HAVE_IMAGEMAGICK == 6)

    case Rec601LumaColorspace:
      return "rec601luma";

    case Rec709LumaColorspace:
      return "rec709luma";

#else

    case LinearGRAYColorspace:
      return "lgray";

    case xyYColorspace:
      return "xyy";

#endif

    case scRGBColorspace:
      return "scrgb";

    case sRGBColorspace:
      return "srgb";

    case RGBColorspace:
      return "rgb" ;

    case TransparentColorspace:
      return "transparent";

    case YCbCrColorspace:
      return "ycbcr";

    case YCCColorspace:
      return "ycc";

    case YDbDrColorspace:
      return "ydbdr";

    case YIQColorspace:
      return "yiq";

    case YPbPrColorspace:
      return "ypbpr";

    case YUVColorspace:
      return "yuv";

    default:
      return "undef";
  }
}

static MetricType string2metric(char * val)
{
  if (!strcmp(val,"avgabs") || !strcmp(val,"meanabs"))
    return MeanAbsoluteErrorMetric;

  if (!strcmp(val,"abs"))
    return AbsoluteErrorMetric;

  if (!strcmp(val,"fuzz"))
    return FuzzErrorMetric;

  if (!strcmp(val,"meansquare") || !strcmp(val,"mse"))
    return MeanSquaredErrorMetric;

  if (!strcmp(val,"rms"))
    return RootMeanSquaredErrorMetric;

  if (!strcmp(val,"peak") || !strcmp(val,"peakabs"))
    return PeakAbsoluteErrorMetric;

  if (!strcmp(val,"normcrosscorr"))
    return NormalizedCrossCorrelationErrorMetric;

#if (HAVE_IMAGEMAGICK > 6)

  if (!strcmp(val,"strsim"))
    return StructuralSimilarityErrorMetric;

  if (!strcmp(val,"strdiss"))
    return StructuralDissimilarityErrorMetric;

  if (!strcmp(val,"peaksnr"))
    return PeakSignalToNoiseRatioErrorMetric;

  if (!strcmp(val,"mepp"))
    return MeanErrorPerPixelErrorMetric;

  if (!strcmp(val,"perchash"))
    return PerceptualHashErrorMetric;

#endif

  return UndefinedErrorMetric;
}

static ChannelType string2channel(char *c)
{
  ChannelType rval=UndefinedChannel;

  if (isvalidstring(c)<1)
    rval |=  DefaultChannels;

  if (strstr(c, "r") != 0)
    rval |=  RedChannel;

  if (strstr(c, "c") != 0)
    rval |=  CyanChannel;

  if (strstr(c, "g") != 0)
    rval |=  GreenChannel;

  if (strstr(c, "m") != 0)
    rval |=  MagentaChannel;

  if (strstr(c, "b") != 0)
    rval |=  BlueChannel;

  if (strstr(c, "y") != 0)
    rval |=  YellowChannel;

  if (strstr(c, "o") != 0)
    rval |=  OpacityChannel;
  if (strstr(c, "a") != 0)
    rval |=  AlphaChannel;

  if (strstr(c, "B") != 0 || strstr(c, "k") != 0)
    rval |=  BlackChannel;

#if (HAVE_IMAGEMAGICK == 6)

  if (strstr(c, "M") != 0 || strstr(c, "t") != 0)
    rval |=  MatteChannel;

#endif

  if (strstr(c, "G") != 0)
    rval |= GrayChannel;

  return rval;
}

CompositeOperator string2compositeoperator(char *c)
{
  if (isvalidstring(c)<1)
    return OverCompositeOp;

  if (!strcmp(c,"mod"))
    return ModulusAddCompositeOp;

  if (!strcmp(c,"blend"))
    return BlendCompositeOp;

  if (!strcmp(c,"changemask"))
    return ChangeMaskCompositeOp;

  if (!strcmp(c,"clear"))
    return ClearCompositeOp;

  if (!strcmp(c,"burn"))
    return ColorBurnCompositeOp;

  if (!strcmp(c,"dodge"))
    return ColorDodgeCompositeOp;

  if (!strcmp(c,"colorize"))
    return ColorizeCompositeOp;

  if (!strcmp(c,"copyblack"))
    return CopyBlackCompositeOp;

  if (!strcmp(c,"copyblue"))
    return CopyBlueCompositeOp;

  if (!strcmp(c,"copycyan"))
    return CopyCyanCompositeOp;

  if (!strcmp(c,"copygreen"))
    return CopyGreenCompositeOp;

  if (!strcmp(c,"copymag"))
    return CopyMagentaCompositeOp;

#if (HAVE_IMAGEMAGICK == 6)

  if (!strcmp(c,"copyopac"))
    return CopyOpacityCompositeOp;

#endif

  if (!strcmp(c,"copyred"))
    return CopyRedCompositeOp;

  if (!strcmp(c,"copyyell"))
    return CopyYellowCompositeOp;

  if (!strcmp(c,"darken"))
    return DarkenCompositeOp;

  if (!strcmp(c,"dstatop"))
    return DstAtopCompositeOp;

  if (!strcmp(c,"dst"))
    return DstCompositeOp;

  if (!strcmp(c,"dstin"))
    return DstInCompositeOp;

  if (!strcmp(c,"dstout"))
    return DstOutCompositeOp;

  if (!strcmp(c,"dstover"))
    return DstOverCompositeOp;

  if (!strcmp(c,"dissolv"))
    return DissolveCompositeOp;

  if (!strcmp(c,"excl"))
    return ExclusionCompositeOp;

  if (!strcmp(c,"hardlight"))
    return HardLightCompositeOp;

  if (!strcmp(c,"huecomp"))
    return HueCompositeOp;

  if (!strcmp(c,"lighten"))
    return LightenCompositeOp;

  if (!strcmp(c,"linlight"))
    return LinearLightCompositeOp;

  if (!strcmp(c,"luminize"))
    return LuminizeCompositeOp;

  if (!strcmp(c,"minusdst"))
    return MinusDstCompositeOp;

  if (!strcmp(c,"modulate"))
    return ModulateCompositeOp;

  if (!strcmp(c,"multiply"))
    return MultiplyCompositeOp;

  if (!strcmp(c,"overlay"))
    return OverlayCompositeOp;

  if (!strcmp(c,"replace"))
    return ReplaceCompositeOp;

  if (!strcmp(c,"saturate"))
    return SaturateCompositeOp;

  if (!strcmp(c,"screen"))
    return ScreenCompositeOp;

  if (!strcmp(c,"softlight"))
    return SoftLightCompositeOp;

  if (!strcmp(c,"srcatop"))
    return SrcAtopCompositeOp;

  if (!strcmp(c,"src"))
    return SrcCompositeOp;

  if (!strcmp(c,"srcin"))
    return SrcInCompositeOp;

  if (!strcmp(c,"srcout"))
    return SrcOutCompositeOp;

  if (!strcmp(c,"srcover"))
    return SrcOverCompositeOp;

  if (!strcmp(c,"modsub"))
    return ModulusSubtractCompositeOp;

  if (!strcmp(c,"threshold"))
    return ThresholdCompositeOp;

  if (!strcmp(c,"divide"))
    return DivideDstCompositeOp;

  if (!strcmp(c,"distort"))
    return DistortCompositeOp;

  if (!strcmp(c,"blur"))
    return BlurCompositeOp;

  if (!strcmp(c,"pegtop"))
    return PegtopLightCompositeOp;

  if (!strcmp(c,"vividlight"))
    return VividLightCompositeOp;

  if (!strcmp(c,"pinlight"))
    return PinLightCompositeOp;

  if (!strcmp(c,"lindodge"))
    return LinearDodgeCompositeOp;

  if (!strcmp(c,"linburn"))
    return LinearBurnCompositeOp;

  if (!strcmp(c,"math"))
    return MathematicsCompositeOp;

  if (!strcmp(c,"divsrc"))
    return DivideSrcCompositeOp;

  if (!strcmp(c,"minsrc"))
    return MinusSrcCompositeOp;

  if (!strcmp(c,"darkenint"))
    return DarkenIntensityCompositeOp;

  if (!strcmp(c,"lightenint"))
    return LightenIntensityCompositeOp;

  if (!strcmp(c,"over"))
    return OverCompositeOp;

  if (!strcmp(c,"in"))
    return InCompositeOp;

  if (!strcmp(c,"out"))
    return OutCompositeOp;

  if (!strcmp(c,"atop"))
    return AtopCompositeOp;

  if (!strcmp(c,"xor"))
    return XorCompositeOp;

  if (!strcmp(c,"plus"))
    return PlusCompositeOp;

  if (!strcmp(c,"diff"))
    return DifferenceCompositeOp;

  if (!strcmp(c,"bumpmap"))
    return BumpmapCompositeOp;

  if (!strcmp(c,"copy"))
    return CopyCompositeOp;

  if (!strcmp(c,"displace"))
    return DisplaceCompositeOp;

  return OverCompositeOp;
}

#if (HAVE_IMAGEMAGICK == 6)
DistortImageMethod
#else
DistortMethod
#endif
               string2distortmethod(char *c)
{
  if (isvalidstring(c)<1)
    return UndefinedDistortion;

  if (!strcmp(c,"affine"))
      return AffineDistortion;

  if (!strcmp(c,"affineproj"))
    return AffineProjectionDistortion;

  if (!strcmp(c,"scalerottran"))
    return ScaleRotateTranslateDistortion;

  if (!strcmp(c,"perspectivedist"))
    return PerspectiveDistortion;

  if (!strcmp(c,"perspectiveproj"))
    return PerspectiveProjectionDistortion;

  if (!strcmp(c,"bilinearforward"))
    return BilinearForwardDistortion;

  if (!strcmp(c,"bilinear"))
    return BilinearDistortion;

  if (!strcmp(c,"bilinearrev"))
    return BilinearReverseDistortion;

  if (!strcmp(c,"poly"))
    return PolynomialDistortion;

  if (!strcmp(c,"arc"))
    return ArcDistortion;

  if (!strcmp(c,"polar"))
    return PolarDistortion;

  if (!strcmp(c,"depolar"))
    return DePolarDistortion;

  if (!strcmp(c,"cyl"))
    return Cylinder2PlaneDistortion;

  if (!strcmp(c,"plane"))
    return Plane2CylinderDistortion;

  if (!strcmp(c,"barrel"))
    return BarrelDistortion;

  if (!strcmp(c,"barrelinv"))
    return BarrelInverseDistortion;

  if (!strcmp(c,"shep"))
    return ShepardsDistortion;

  if (!strcmp(c,"resize"))
    return ResizeDistortion;

  if (!strcmp(c,"sent"))
    return SentinelDistortion ;

  return UndefinedDistortion;
}

#if (HAVE_IMAGEMAGICK == 6)
InterpolatePixelMethod
#else
PixelInterpolateMethod
#endif
    string2pixelinterpolatemethod(char *c)
{
  if (isvalidstring(c)<1)
    return UndefinedInterpolatePixel;

  if (!strcmp(c,"avg"))
    return AverageInterpolatePixel;

  /* Triangular/bilinear filter interpolation */
  if (!strcmp(c,"bilin") || !strcmp(c,"tri"))
    return BilinearInterpolatePixel;

  /* Integer (floor) interpolation */
  if (!strcmp(c,"int") || !strcmp(c,"floor"))
    return IntegerInterpolatePixel;

  /* Cubic Spline (blurred) interpolation */
  if (!strcmp(c,"spl"))
    return SplineInterpolatePixel;

#if (HAVE_IMAGEMAGICK == 6)

  if (!strcmp(c,"bicubic"))
    return BicubicInterpolatePixel;

  if (!strcmp(c,"filter"))
    return FilterInterpolatePixel;

#else

  /* Nearest Neighbour Only */
  if (!strcmp(c,"nn"))
    return NearestInterpolatePixel;

  /* Average 4 nearest neighbours */
  if (!strcmp(c,"avg:4"))
    return AverageInterpolatePixel;

  /* Average 9 nearest neighbours */
  if (!strcmp(c,"avg:9"))
    return Average9InterpolatePixel;

  /* Average 16 nearest neighbours */
  if (!strcmp(c,"avg:16"))
    return Average16InterpolatePixel;

  /* Just return background color */
  if (!strcmp(c,"bkg"))
    return BackgroundInterpolatePixel;

  /* blend of nearest 1, 2 or 4 pixels */
  if (!strcmp(c,"blend"))
    return BlendInterpolatePixel;

  /* Catmull-Rom interpolation */
  if (!strcmp(c,"catrom"))
    return CatromInterpolatePixel;

  /* Triangular Mesh interpolation */
  if (!strcmp(c,"mesh"))
    return MeshInterpolatePixel;

#endif

  return UndefinedInterpolatePixel;
}


static const char * interpolatepixelmethod2string(
#if (HAVE_IMAGEMAGICK == 6)
    InterpolatePixelMethod
#else
    PixelInterpolateMethod
#endif
        interpolate)
{
  switch (interpolate)
  {

    /* Triangular/bilinear filter interpolation */
    case BilinearInterpolatePixel:
      return "bilin";

    /* Integer (floor) interpolation */
    case IntegerInterpolatePixel:
      return "int";

    /* Cubic Spline (blurred) interpolation */
    case SplineInterpolatePixel:
      return "spl";

#if (HAVE_IMAGEMAGICK == 6)

    case AverageInterpolatePixel:
      return "avg";

    case BicubicInterpolatePixel:
      return "bicubic";

    case FilterInterpolatePixel:
      return "filter";

#else

    /* Nearest Neighbour Only */
    case NearestInterpolatePixel:
      return "nn";

    /* Average 4 nearest neighbours */
    case  AverageInterpolatePixel:
      return "avg:4";

    /* Average 9 nearest neighbours */
    case  Average9InterpolatePixel:
      return "avg:9";

    /* Average 16 nearest neighbours */
    case  Average16InterpolatePixel:
      return "avg:16";

    /* Just case  background color */
    case  BackgroundInterpolatePixel:
      return "bkg";

    /* blend of nearest 1, 2 or 4 pixels */
    case  BlendInterpolatePixel:
      return "blend";

    /* Catmull-Rom interpolation */
    case  CatromInterpolatePixel:
      return "catrom";

    /* Triangular Mesh interpolation */
    case  MeshInterpolatePixel:
      return "mesh";

#endif

    default:
      return "undef";
  }
}



MagickEvaluateOperator string2op(char *expr)
{
  if (isvalidstring(expr)<1)
    return UndefinedEvaluateOperator;

  if (!strcmp(expr,"+"))
    return AddEvaluateOperator;

  if (!strcmp(expr,"and"))
    return AndEvaluateOperator;

  if (!strcmp(expr,"/") || !strcmp(expr,"divide"))
    return DivideEvaluateOperator;

  if (!strcmp(expr,"<<") || !strcmp(expr,"shiftleft"))
    return LeftShiftEvaluateOperator;

  if (!strcmp(expr,"max"))
    return MaxEvaluateOperator;

  if (!strcmp(expr,"min"))
    return MinEvaluateOperator;

  if (!strcmp(expr,"*") || !strcmp(expr,"multiply"))
    return MultiplyEvaluateOperator;

  if (!strcmp(expr,"or"))
    return OrEvaluateOperator;

  if (!strcmp(expr,">>") || !strcmp(expr,"shiftright"))
    return RightShiftEvaluateOperator;

  if (!strcmp(expr,"set"))
    return SetEvaluateOperator;

  if (!strcmp(expr,"-") || !strcmp(expr,"subtract") || !strcmp(expr,"minus"))
    return SubtractEvaluateOperator;

  if (!strcmp(expr,"xor"))
    return XorEvaluateOperator;

  if (!strcmp(expr,"pow"))
    return PowEvaluateOperator;

  if (!strcmp(expr,"log"))
    return LogEvaluateOperator;

  if (!strcmp(expr,"threshold"))
    return ThresholdEvaluateOperator;

  if (!strcmp(expr,"thresholdblack"))
    return ThresholdBlackEvaluateOperator;

  if (!strcmp(expr,"thresholdwhite"))
    return ThresholdWhiteEvaluateOperator;

  if (!strcmp(expr,"gaussian") || !strcmp(expr,"noisegaussian"))
    return GaussianNoiseEvaluateOperator;

  if (!strcmp(expr,"impulse") || !strcmp(expr,"noiseimpulse"))
    return ImpulseNoiseEvaluateOperator;

  if (!strcmp(expr,"laplace") || !strcmp(expr,"noiselaplace"))
    return LaplacianNoiseEvaluateOperator;

  if (!strcmp(expr,"multiplicative") || !strcmp(expr,"noisemultiplicative"))
    return MultiplicativeNoiseEvaluateOperator;

  if (!strcmp(expr,"poisson") || !strcmp(expr,"noisepoisson"))
    return PoissonNoiseEvaluateOperator;

  if (!strcmp(expr,"uniform") || !strcmp(expr,"noiseuniform"))
    return UniformNoiseEvaluateOperator;

  if (!strcmp(expr,"cos"))
    return CosineEvaluateOperator;

  if (!strcmp(expr,"sin"))
    return SineEvaluateOperator;

  if (!strcmp(expr,"+mod") || !strcmp(expr,"addmodulus"))
    return AddModulusEvaluateOperator;

  if (!strcmp(expr,"mean"))
    return MeanEvaluateOperator;

  if (!strcmp(expr,"abs"))
    return AbsEvaluateOperator;

  if (!strcmp(expr,"exp"))
    return ExponentialEvaluateOperator;

  if (!strcmp(expr,"median"))
    return MedianEvaluateOperator;

  return UndefinedEvaluateOperator;
}

MagickFunction string2func(char *expr)
{
  if (isvalidstring(expr)<1)
    return UndefinedFunction;

  if (!strcmp(expr,"asin") || !strcmp(expr,"arcsin"))
    return ArcsinFunction;

  if (!strcmp(expr,"atan") || !strcmp(expr,"arctan"))
    return ArctanFunction;

  if (!strcmp(expr,"poly"))
    return PolynomialFunction;

  if (!strcmp(expr,"sin"))
    return SinusoidFunction;

  return UndefinedFunction;
}

InterlaceType string2interlace(char *val)
{
  if (isvalidstring(val) < 1)
    return NoInterlace;

  if (!strcmp(val,"line"))
    return LineInterlace;

   if (!strcmp(val,"plane"))
    return PlaneInterlace;

   if (!strcmp(val,"part"))
    return PartitionInterlace;

   if (!strcmp(val,"gif"))
    return GIFInterlace;

   if (!strcmp(val,"jpeg"))
    return JPEGInterlace;

   if (!strcmp(val,"png"))
    return PNGInterlace;

  return NoInterlace;
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
  if (columns < 1)
    c = (size_t) ((double) MagickGetImageWidth(wand) * columns);
  else
    c = (size_t) columns;
  if (rows <= 1)
    r = (size_t) ((double) MagickGetImageHeight(wand) * rows);
  else
    r = (size_t) rows;
  return MagickAdaptiveResizeImage((MagickWand *)wand, c, r);
}

static MagickBooleanType
    _MagickAdaptiveThresholdImage(MagickWand *wand, double width,double height,double offset)
{
  return MagickAdaptiveThresholdImage((MagickWand *)wand,
                                       (const size_t) width,
                                        (const size_t) height,
                                         (const size_t) offset);
}

static MagickBooleanType
    _MagickCycleColormapImage(MagickWand *wand, double displace)
{
  return MagickCycleColormapImage((MagickWand *)wand, (const ssize_t) displace);
}

static MagickBooleanType
    _MagickContrastImage(MagickWand *wand, double sharpen)
{
  return MagickContrastImage((MagickWand *)wand, (const MagickBooleanType) sharpen);
}

static MagickBooleanType
    _MagickNegateImage(MagickWand *wand, double gray)
{
  return MagickNegateImage((MagickWand *)wand, (const MagickBooleanType) gray);
};

static MagickBooleanType
    _MagickPosterizeImage(MagickWand *wand,double levels,double dither)
{
  return MagickPosterizeImage((MagickWand *)wand, (const size_t) levels,
                               (const MagickBooleanType) dither);
}

static MagickBooleanType
    _MagickRollImage(MagickWand *wand, double x, double y)
{
  return MagickRollImage((MagickWand *)wand, (const size_t) x, (const size_t) y);
}

static MagickBooleanType
    _MagickSampleImage(MagickWand *wand, double columns, double rows)
{
  return MagickSampleImage((MagickWand *)wand, (const size_t) columns, (const size_t) rows);
}

static MagickBooleanType
    _MagickScaleImage(MagickWand *wand, double columns, double rows)
{
  return MagickScaleImage((MagickWand *)wand, (const size_t) columns, (const size_t) rows);
}

static MagickBooleanType
    _MagickShadeImage(MagickWand *wand, double gray, double azimuth, double elevation)
{
  return MagickShadeImage((MagickWand *)wand, (const MagickBooleanType) gray,
                           (const double) azimuth, (const double) elevation);
}

static MagickBooleanType
    _MagickShadowImage(MagickWand *wand, double opacity, double sigma, double x, double y)
{
  return MagickShadowImage((MagickWand *)wand, (const double) opacity,
                            (const double)sigma, (const size_t) x,
                             (const size_t) y);
}

static MagickBooleanType
    _MagickShaveImage(MagickWand *wand, double columns, double rows)
{
  return MagickShaveImage((MagickWand *)wand, (const size_t) columns, (const size_t) rows);
}

static MagickBooleanType
    _MagickSigmoidalContrastImage(MagickWand *wand, double sharpen, double alpha, double beta)
{
  return MagickSigmoidalContrastImage((MagickWand *) wand,
                                       (const MagickBooleanType) sharpen,
                                        (const double) alpha,
                                         (const double) beta);
}

static MagickBooleanType
    _MagickSpliceImage(MagickWand *wand, double width, double height, double x, double y)
{
  return MagickSpliceImage((MagickWand *)wand,
                            (const size_t) width, (const size_t) height, (const size_t) x,
                             (const size_t) y);
}

static MagickBooleanType
    _MagickThumbnailImage(MagickWand *wand, double columns, double rows)
{
  return MagickThumbnailImage((MagickWand *)wand,(const size_t) columns,(const size_t) rows);
}

static MagickBooleanType
    _MagickVignetteImage(MagickWand *wand, double radius, double sigma, double x,double y)
{
  return MagickVignetteImage((MagickWand *)wand,
                              (const double) radius,
                               (const double) sigma,
                                (const size_t) x,
                                 (const size_t) y);
}


#if (HAVE_IMAGEMAGICK == 6)
#include "r_imagemagick6.c"
#endif

#if (HAVE_IMAGEMAGICK == 7)
#include "r_imagemagick7.c"
#endif

//
// set default wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".wand"
Ent * ent_im_wand (int nargs, Datum args[])
{
  Ent *e1=0, *rent, *R;
  int i,k,icurr;
  short int *image_idx = (short int *) string_buff;
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


#undef  THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".write"
Ent * ent_im_wand_write_image (int nargs, Datum args[])
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
Ent * ent_im_artifacts (int nargs, Datum args[])
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
      sprintf(string_buff,"%g", dummy);
      new_val = string_buff;
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
Ent * ent_im_options (int nargs, Datum args[])
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
      sprintf(string_buff,"%g", dummy);
      new_val = string_buff;
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
#include <sys/wait.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <linux/sched.h>
#include <sched.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#define _XOPEN_SOURCE 700
#include <unistd.h>
static int pid3=0;

#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".disp"
Ent * ent_im_display_wand (int nargs, Datum args[])
{
  Ent *e1=0;
  char *c=0, *xserv;
  int status;

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
  if (isvalidstring(c)>0)
  {
    xserv = c;
  }
  else
  {
    xserv = getenv("DISPLAY");
  }

  pid_t p1 = fork();
  if (p1 != 0)
  {
    wait(&status);
  }
  else
  {
    pid_t p2 = fork();
    if (p2 == 0)
    {
      pid3 = getpid();
      MagickDisplayImage(magick_wand[default_wand_idx], xserv);
    }

    exit(0);
  }


  ent_Clean(e1);
  return ent_Create_Rlab_Success();
}

Ent * ent_im_display_wand_working (int nargs, Datum args[])
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
  if (isvalidstring(c)>0)
  {
    xsrv = c;
  }
  else
  {
    xsrv = getenv("DISPLAY");
  }

  MagickDisplayImage(magick_wand[default_wand_idx], xsrv);

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
Ent * ent_im_read_image_to_wand (int nargs, Datum args[])
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
Ent * ent_im_wand_iterate_images (int nargs, Datum args[])
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
Ent * ent_im_add (int nargs, Datum args[])
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
Ent * ent_im_clone (int nargs, Datum args[])
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
Ent * ent_im_property (int nargs, Datum args[])
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
Ent * ent_im_clear (int nargs, Datum args[])
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
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".close"
Ent * ent_im_close (int nargs, Datum args[])
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
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".exit"
Ent * ent_im_exit (int nargs, Datum args[])
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
Ent * ent_im_wand_join (int nargs, Datum args[])
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
             THIS_SOLVER ":   'another_wand' contains images that are inserted into default "
                 "Magick Wand depending on the iterator.\n");
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
Ent * ent_im_wand_append (int nargs, Datum args[])
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
             THIS_SOLVER ": Append images in default Magick Wand from active image onwards to "
                 "a new Magick Wand,\n");
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
// deconstruct default wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".deconstruct"
Ent * ent_im_wand_deconstruct (int nargs, Datum args[])
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
Ent * ent_im_wand_distort (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *x=0;
  char *c=0;
  size_t n=0;

#if (HAVE_IMAGEMAGICK == 6)
  DistortImageMethod
#else
  DistortMethod
#endif
      method=AffineDistortion;

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
    method = string2distortmethod(c);
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
#if (HAVE_IMAGEMAGICK == 6)
          PixelSetOpacity(bkgrnd, qrange);
#else
          PixelSetAlpha(bkgrnd, qrange);
#endif

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
              pixels_16[i*npf*width + npf*j + k] =
                  (unsigned short int) Mdi0(pixel_data[k],i,j) > (unsigned short int) qrange
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


#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".specs"
Ent * ent_im_wand_image_specs(int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  char *spec=0;
  int i;
  double d;
  MDR *w=0;

  if (!magick_wand[default_wand_idx] || !InitMagickWandLib)
    rerror(THIS_SOLVER  ": Help! Library not initialized!\n");

  if ((nargs != 0) && (nargs != 2))
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
      if (d<-1 || d>1)
        rerror(THIS_SOLVER ": Second argument 'val' must be integer or real scalar!");

      if (d==-1)
      {
        ent_data(rent) = mdi_CreateScalar( MagickSetImageAlphaChannel(magick_wand[default_wand_idx],
                                  DeactivateAlphaChannel) );
      }
      else
      {
        ent_data(rent) = mdi_CreateScalar( MagickSetImageAlphaChannel(magick_wand[default_wand_idx],
                 ActivateAlphaChannel) );
        if (d == 0)
        {
          ent_data(rent) = mdi_CreateScalar( MagickSetImageAlphaChannel(magick_wand[default_wand_idx],
                   OpaqueAlphaChannel) );
        }
        else if (d == 1)
        {
          ent_data(rent) = mdi_CreateScalar( MagickSetImageAlphaChannel(magick_wand[default_wand_idx],
                   SetAlphaChannel) );
        }
        else
        {
#if (HAVE_IMAGEMAGICK == 6)
          ent_data(rent) = mdi_CreateScalar( MagickSetImageOpacity(magick_wand[default_wand_idx], d) );
#else
          ent_data(rent) = mdi_CreateScalar( MagickSetImageAlpha(magick_wand[default_wand_idx], d) );
#endif
        }
      }
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
      if (SIZE(w) != 4)
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
      if (SIZE(w) < 3)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      PixelWand *pw = NewPixelWand();
      PixelSetRed     (pw, mdrV0(w,0));
      PixelSetGreen   (pw, mdrV0(w,1));
      PixelSetBlue    (pw, mdrV0(w,2));
#if (HAVE_IMAGEMAGICK == 6)
      if (SIZE(w) == 4)
        PixelSetOpacity (pw, MdrV0(w,3));
      else
        PixelSetOpacity (pw, 0.0);
#else
      if (SIZE(w) == 4)
        PixelSetAlpha (pw, mdrV0(w,3));
      else
        PixelSetAlpha (pw, 1.0);
#endif
      ent_data(rent) = mdi_CreateScalar( MagickSetImageBackgroundColor(magick_wand[default_wand_idx],
               pw) );
      pw = DestroyPixelWand(pw);
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
#if (HAVE_IMAGEMAGICK == 6)
    else if (!strcmp(spec, "bias"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real scalar!");
      d = (int) class_double(e2);

      ent_data(rent) = mdi_CreateScalar( MagickSetImageBias(magick_wand[default_wand_idx], d) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
#endif
    else if (!strcmp(spec, "blue_primary"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
#if (HAVE_IMAGEMAGICK == 6)
      if (SIZE(w)!=2)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      ent_data(rent) = mdi_CreateScalar( MagickSetImageBluePrimary(magick_wand[default_wand_idx],
               mdrV0(w,0), mdrV0(w,1)) );
#else
      if (SIZE(w)!=3)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      ent_data(rent) = mdi_CreateScalar( MagickSetImageBluePrimary(magick_wand[default_wand_idx],
               mdrV0(w,0), mdrV0(w,1), mdrV0(w,2)) );
#endif
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "red_primary"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
#if (HAVE_IMAGEMAGICK == 6)
      if (SIZE(w)!=2)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      ent_data(rent) = mdi_CreateScalar( MagickSetImageRedPrimary(magick_wand[default_wand_idx],
               mdrV0(w,0), mdrV0(w,1)) );
#else
      if (SIZE(w)!=3)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      ent_data(rent) = mdi_CreateScalar( MagickSetImageRedPrimary(magick_wand[default_wand_idx],
               mdrV0(w,0), mdrV0(w,1), mdrV0(w,2)) );
#endif
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "green_primary"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
#if (HAVE_IMAGEMAGICK == 6)
      if (SIZE(w)!=2)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      ent_data(rent) = mdi_CreateScalar( MagickSetImageGreenPrimary(magick_wand[default_wand_idx],
               mdrV0(w,0), mdrV0(w,1)) );
#else
      if (SIZE(w)!=3)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      ent_data(rent) = mdi_CreateScalar( MagickSetImageGreenPrimary(magick_wand[default_wand_idx],
               mdrV0(w,0), mdrV0(w,1), mdrV0(w,2)) );
#endif
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "border_color"))
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      w = class_matrix_real(e2);
      if (SIZE(w)<3)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");

      PixelWand *pw = NewPixelWand();
      PixelSetRed     (pw, (double) MdrV0(w,0));
      PixelSetGreen   (pw, (double) MdrV0(w,1));
      PixelSetBlue    (pw, (double) MdrV0(w,2));
#if (HAVE_IMAGEMAGICK == 6)
      if (SIZE(w) == 4)
        PixelSetOpacity (pw, MdrV0(w,3));
      else
        PixelSetOpacity (pw, 0.0);
#else
      if (SIZE(w) == 4)
        PixelSetAlpha (pw, mdrV0(w,3));
      else
        PixelSetAlpha (pw, 1.0);
#endif
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
      if (SIZE(w)<3)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");

      PixelWand *pw = NewPixelWand();
      PixelSetRed     (pw, (double) MdrV0(w,0));
      PixelSetGreen   (pw, (double) MdrV0(w,1));
      PixelSetBlue    (pw, (double) MdrV0(w,2));
#if (HAVE_IMAGEMAGICK == 6)
      if (SIZE(w) == 4)
        PixelSetOpacity (pw, MdrV0(w,3));
      else
        PixelSetOpacity (pw, 0.0);
#else
      if (SIZE(w) == 4)
        PixelSetAlpha (pw, mdrV0(w,3));
      else
        PixelSetAlpha (pw, 1.0);
#endif
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
      cs = string2colorspace(val);
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

      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");
      val = class_char_pointer(e2);
      InterlaceType interlace=string2interlace(val);
      ent_data(rent) = mdi_CreateScalar( MagickSetImageInterlaceScheme( magick_wand[default_wand_idx],
               interlace) );
      ent_type(rent) = MATRIX_DENSE_REAL;
    }
    else if (!strcmp(spec, "interpolate"))
    {
      // Oh gravity! It's pulling me!
      char *val=0;
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");
      val = class_char_pointer(e2);
#if (HAVE_IMAGEMAGICK == 6)
      InterpolatePixelMethod
#else
      PixelInterpolateMethod
#endif
          interpolate=string2pixelinterpolatemethod(val);
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
      if ((SIZE(w) != 1) && (SIZE(w) != 3) && (SIZE(w) != 4))
        rerror(THIS_SOLVER ": Second argument 'val' must be integer or real scalar!");

      if (SIZE(w) == 1)
      {
        if (mdrV0(w,0))
          matte = MagickTrue;
        ent_data(rent) = mdi_CreateScalar( MagickSetImageMatte(magick_wand[default_wand_idx],matte) );
      }
      else
      {
        PixelWand *pw = NewPixelWand();
        PixelSetRed     (pw, mdrV0(w,0));
        PixelSetGreen   (pw, mdrV0(w,1));
        PixelSetBlue    (pw, mdrV0(w,2));
#if (HAVE_IMAGEMAGICK == 6)
        if (SIZE(w) == 4)
          PixelSetOpacity (pw, mdrV0(w,3));
        else
          PixelSetOpacity (pw, 0.0);
#else
        if (SIZE(w) == 4)
          PixelSetAlpha (pw, mdrV0(w,3));
        else
          PixelSetAlpha (pw, 1.0);
#endif
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

      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) != MATRIX_DENSE_STRING)
        rerror(THIS_SOLVER ": Second argument 'val' must be string!");
      val = class_char_pointer(e2);
      ImageType imagetype=string2imagetype(val);
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
#if (HAVE_IMAGEMAGICK == 6)
      if (SIZE(w)!=2)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      ent_data(rent) = mdi_CreateScalar( MagickSetImageWhitePoint(magick_wand[default_wand_idx],
               mdrV0(w,0), mdrV0(w,1)) );
#else
      if (SIZE(w)!=3)
        rerror(THIS_SOLVER ": Second argument 'val' must be real vector!");
      ent_data(rent) = mdi_CreateScalar( MagickSetImageWhitePoint(magick_wand[default_wand_idx],
               mdrV0(w,0), mdrV0(w,1), mdrV0(w,2)) );
#endif
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
#if (HAVE_IMAGEMAGICK == 6)
  MDR *r = mdr_Create(1,2);
  MagickGetImageRedPrimary(magick_wand[default_wand_idx], &MdrV0(r,0), &MdrV0(r,1));
#else
  MDR *r = mdr_Create(1,3);
  MagickGetImageRedPrimary(magick_wand[default_wand_idx], &MdrV0(r,0), &MdrV0(r,1), &MdrV0(r,2));
#endif
  E = ent_Create ();
  ent_data (E) = r;
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("red_primary"), E);

  // green primary:
#if (HAVE_IMAGEMAGICK == 6)
  MDR *g = mdr_Create(1,2);
  MagickGetImageGreenPrimary(magick_wand[default_wand_idx], &MdrV0(r,0), &MdrV0(r,1));
#else
  MDR *g = mdr_Create(1,3);
  MagickGetImageGreenPrimary(magick_wand[default_wand_idx], &MdrV0(r,0), &MdrV0(r,1), &MdrV0(r,2));
#endif
  E = ent_Create ();
  ent_data (E) = g;
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("green_primary"), E);

  // blue primary:
#if (HAVE_IMAGEMAGICK == 6)
  MDR *b = mdr_Create(1,2);
  MagickGetImageBluePrimary(magick_wand[default_wand_idx], &MdrV0(r,0), &MdrV0(r,1));
#else
  MDR *b = mdr_Create(1,3);
  MagickGetImageBluePrimary(magick_wand[default_wand_idx], &MdrV0(r,0), &MdrV0(r,1), &MdrV0(r,2));
#endif
  E = ent_Create ();
  ent_data (E) = b;
  ent_type (E) = MATRIX_DENSE_REAL;
  install (bw, ("blue_primary"), E);

  // colorspace
  ColorspaceType cs = MagickGetImageColorspace( magick_wand[default_wand_idx] );
  E = ent_Create ();
  MDS * sp = mds_CreateScalar((char *) colorspace2string(cs));
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
#if (HAVE_IMAGEMAGICK == 6)
  InterpolatePixelMethod
#else
  PixelInterpolateMethod
#endif
      interpolate = MagickGetImageInterpolateMethod(magick_wand[default_wand_idx]);
  MDS *interp=mds_CreateScalar((char *)interpolatepixelmethod2string(interpolate));
  RLABCODE_INSTALL_ENTITY_IN_BTREE(bw,E,MATRIX_DENSE_STRING,interp,"interpolate");

  // orientation
#if (HAVE_IMAGEMAGICK == 6)
  InterpolatePixelMethod
#else
  PixelInterpolateMethod
#endif
      orientation = MagickGetImageInterpolateMethod(magick_wand[default_wand_idx]);
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
  RLABCODE_INSTALL_ENTITY_IN_BTREE(bw,E,MATRIX_DENSE_STRING,orient,"orientation");

  // page: on 64-bit machines size(size_t) > size(int)
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
  RLABCODE_INSTALL_ENTITY_IN_BTREE(bw,E,MATRIX_DENSE_REAL,pag,"page");

  // resolution
  MDR *resol  = mdr_Create(1,2);
  MagickGetImageResolution(magick_wand[default_wand_idx], &MdrV0(resol,0), &MdrV0(resol,1));
  RLABCODE_INSTALL_ENTITY_IN_BTREE(bw,E,MATRIX_DENSE_REAL,resol,"resolution");

  // for what follows:
  PixelWand *pw = NewPixelWand();
  MDR *rgba=0;

  // background
  MagickGetImageBackgroundColor(magick_wand[default_wand_idx], pw);
  rgba = mdr_Create(1,4);
  MdrV0(rgba,0) = PixelGetRed     (pw);
  MdrV0(rgba,1) = PixelGetGreen   (pw);
  MdrV0(rgba,2) = PixelGetBlue    (pw);
#if (HAVE_IMAGEMAGICK == 6)
  MdrV0(rgba,3) = PixelGetOpacity (pw);
#else
  MdrV0(rgba,3) = PixelGetAlpha(pw);
#endif
  RLABCODE_INSTALL_ENTITY_IN_BTREE(bw,E,MATRIX_DENSE_REAL,rgba,"background");

  // border color
  MagickGetImageBorderColor(magick_wand[default_wand_idx], pw);
  rgba = mdr_Create(1,4);
  MdrV0(rgba,0) = PixelGetRed     (pw);
  MdrV0(rgba,1) = PixelGetGreen   (pw);
  MdrV0(rgba,2) = PixelGetBlue    (pw);
#if (HAVE_IMAGEMAGICK == 6)
  MdrV0(rgba,3) = PixelGetOpacity (pw);
#else
  MdrV0(rgba,3) = PixelGetAlpha(pw);
#endif
  RLABCODE_INSTALL_ENTITY_IN_BTREE(bw,E,MATRIX_DENSE_REAL,rgba,"border_color");

  pw = DestroyPixelWand(pw);

  // image type
  ImageType imval=MagickGetImageType(magick_wand[default_wand_idx]);
  MDS *imtype=mds_CreateScalar((char *) imagetype2string(imval));
  RLABCODE_INSTALL_ENTITY_IN_BTREE(bw,E,MATRIX_DENSE_STRING,imtype,"image_type");

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
  RLABCODE_INSTALL_ENTITY_IN_BTREE(bw,E,MATRIX_DENSE_STRING,us,"resolution_unit");

  // get image white point
#if (HAVE_IMAGEMAGICK == 6)
  MDR *whp = mdr_Create(1,2);
  MagickGetImageWhitePoint(magick_wand[default_wand_idx], &MdrV0(whp,0),&MdrV0(whp,1));
#else
  MDR *whp = mdr_Create(1,3);
  MagickGetImageWhitePoint(magick_wand[default_wand_idx], &MdrV0(whp,0),&MdrV0(whp,1),&MdrV0(whp,2));
#endif
  RLABCODE_INSTALL_ENTITY_IN_BTREE(bw,E,MATRIX_DENSE_REAL,whp,"white_point");

  rent = ent_Create ();
  ent_data (rent) = bw;
  ent_type (rent) = BTREE;
  return rent;
}



#endif  /* HAVE_IMAGEMAGICK */
