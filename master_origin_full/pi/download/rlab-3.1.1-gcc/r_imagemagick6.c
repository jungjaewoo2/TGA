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

#define THIS_FILE "r_imagemagick6.c"

//
// wrapture is here ! repent!
//    Or maybe not. Just wrappers for the functions so that they accept double arguments:
//
static MagickBooleanType
    _MagickNegateImageChannel(MagickWand *wand, const ChannelType channel, double gray)
{
  return MagickNegateImageChannel((MagickWand *)wand, (const ChannelType) channel,
                                  (const MagickBooleanType) gray);
}

static MagickBooleanType
    _MagickPolaroidImage(MagickWand *wand,double angle)
{
  MagickBooleanType rval;
  DrawingWand *drawing_wand = NewDrawingWand();
  rval = MagickPolaroidImage((MagickWand *)wand,
                             (const DrawingWand *) drawing_wand,
                             (const double) angle);
  drawing_wand = DestroyDrawingWand(drawing_wand);
  return rval;
}

static MagickBooleanType
    _MagickSigmoidalContrastImageChannel(MagickWand *wand, ChannelType channel,
                                         double sharpen, double alpha, double beta)
{
  return MagickSigmoidalContrastImageChannel((MagickWand *) wand, (const ChannelType) channel,
                                      (const MagickBooleanType) sharpen,
                                      (const double) alpha,
                                      (const double) beta);
}




//
// set default wand
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".method"
#define IM_MAX_NO_METHOD_PARAMS 12
Ent * ent_im_wand_apply (int nargs, Datum args[])
{
  Ent *e1=0, *en=0, *rent;
  int i,rval=0;
  char *method=0;
  int np=0, have_ch=0;  // for counting parameters to the methods
  char *lb=0,*rb=0, *c=0;
  double params[IM_MAX_NO_METHOD_PARAMS] = {0};
  MDR *mdr_param=0;
  ChannelType ch=UndefinedChannel;

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
    RLABCODE_PROCESS_ARG1_S(THIS_SOLVER,0,e1,method,1);

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
    func_ch = &MagickConvolveImageChannel;
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
      ch = string2channel(lb);
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
      sprintf(string_buff, "%g", class_double(en));
      ent_Clean(en);

      MagickSetImageArtifact(magick_wand[default_wand_idx], "attenuate", string_buff);
    }

    if (ch != UndefinedChannel)
    {
      func_ch(magick_wand[default_wand_idx], ch, noise_type);
    }
    else
    {
      func(magick_wand[default_wand_idx], noise_type);
    }
  }
  else if (func == &MagickResizeImage)
  {
    // special methods:
    if(nargs != 4)
      rerror (THIS_SOLVER  ": " RLAB_ERROR_ARG_4 "\n");

    // second argument is a 1x2 matrix
    RLABCODE_PROCESS_ARG_F_S(THIS_SOLVER,1,en,
                             MATRIX_DENSE_REAL,mdr_param,SIZE,!=,2,
                             RLAB_ERROR_ARG2_VEC_LEN_2
                            );
    params[0] = mdrV0(mdr_param,0);
    params[1] = mdrV0(mdr_param,1);
    ent_Clean(en);

    // third parameter is the filter name
    RLABCODE_PROCESS_ARG_S(THIS_SOLVER,2,en,c,1,RLAB_ERROR_ARG3_MDS_SCALAR);
    FilterTypes filter;
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
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb for the background as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (SIZE(mdr_param) !=3 && SIZE(mdr_param) !=4)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires rgb[a] for"
          " the background as second argument!\n", method);
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
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'angle' as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[4] = class_double(en);
    ent_Clean(en);

    PixelWand *pw = NewPixelWand();
    PixelSetRed(pw, params[0]);
    PixelSetGreen(pw, params[1]);
    PixelSetBlue(pw, params[2]);
    PixelSetOpacity(pw, params[3]);
    MagickSetImageBorderColor(magick_wand[default_wand_idx], pw);
    MagickRotateImage(magick_wand[default_wand_idx], pw, params[4]);
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
    if (SIZE(mdr_param)==4 )
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
    if((nargs != 3) && (nargs != 4))
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
    if (SIZE(mdr_param) !=2)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires [x,y] real vector as third argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    params[0] = mdrV0(mdr_param,0);
    params[1] = mdrV0(mdr_param,1);
    ent_Clean(en); en = 0;

    // third argument:
    //  list with different text properties
    char *font_name=0;
    double angle=0;
    double font_size=0;
    char *font_color=0;
    int aa_flag=1;
    if (nargs == 4)
    {
      en = bltin_get_ent (args[3]);
      if (ent_type(en) == BTREE)
      {
        ListNode *node=0;
        RLABCODE_PROCESS_BTREE_ENTRY_S(en,node,RLAB_NAME_MAGICKWAND_FONT,font_name,6,NULL);
        RLABCODE_PROCESS_BTREE_ENTRY_S(en,node,RLAB_NAME_MAGICKWAND_FONT_COLOR,font_color,1,NULL);
        RLABCODE_PROCESS_BTREE_ENTRY_DOUBLE(en,node,RLAB_NAME_MAGICKWAND_ANGLE,angle);
        RLABCODE_PROCESS_BTREE_ENTRY_DOUBLE(en,node,RLAB_NAME_MAGICKWAND_FONT_SIZE,font_size);
        RLABCODE_PROCESS_BTREE_ENTRY_BOOL(en,node,RLAB_NAME_MAGICKWAND_FONT_ANTIALIAS,aa_flag);
      }
    }
    params[2] = angle;

    DrawingWand *drawing_wand = NewDrawingWand();
    PixelWand   *pixel_wand = NewPixelWand();

    if (font_color)
    {
      PixelSetColor(pixel_wand, font_color);
      DrawSetFillColor(drawing_wand,pixel_wand);
    }
    if (font_name)
      DrawSetFont (drawing_wand, font_name );
    if (font_size)
      DrawSetFontSize(drawing_wand,font_size);
    if (aa_flag)
      DrawSetTextAntialias(drawing_wand,MagickTrue);

    rval = MagickAnnotateImage(magick_wand[default_wand_idx],drawing_wand,
                               params[0],params[1],params[2],ann);

    drawing_wand = DestroyDrawingWand(drawing_wand);
    pixel_wand = DestroyPixelWand(pixel_wand);

    ent_Clean(es);
    ent_Clean(en);
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
      while( ( string_buff[i] = fgetc(fr) ) != EOF  && i<MAX_STRING_BUFF)
        i++;
      if (i==MAX_STRING_BUFF)
      {
        fprintf(stderr, THIS_SOLVER  ": file %s is too big! It should contain a single profile!\n", filename);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      string_buff[i]='\0';
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
      if (SIZE(mdr_param)!=3)
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
      if (SIZE(mdr_param)!=3)
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
        fprintf(stderr, THIS_SOLVER  ": method %s requires 'saturation' real scalar as "
            "fifth argument!\n", method);
        rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
      }
      params[9] = class_double(en);
      ent_Clean(en);

      // now assemble XML file
      i=0;
      char *current_text=string_buff;
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
      fprintf(stderr, "%s", string_buff);

    bval = MagickColorDecisionListImage(magick_wand[default_wand_idx], string_buff);
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
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'color' real"
          " 3- or 4-vector as second argument!\n", method);
    rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
  }
  mdr_param = ent_data(en);
  if ((SIZE(mdr_param)!=3) && (SIZE(mdr_param)!=4))
  {
    fprintf(stderr, THIS_SOLVER  ": method %s requires 'color' real 3- or 4-vector"
        " as second argument!\n", method);
    rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
  }

  PixelWand *c = NewPixelWand();
  PixelSetRed   (c, mdrV0(mdr_param,0));
  PixelSetGreen (c, mdrV0(mdr_param,1));
  PixelSetBlue  (c, mdrV0(mdr_param,2));

  PixelWand *o = NewPixelWand();
  if (SIZE(mdr_param) == 4)
    PixelSetOpacity(o, mdrV0(mdr_param,3));
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
    if (ent_type(en) != MATRIX_DENSE_STRING)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires string scalar describing real square matrix"
          " dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    char * ks = class_char_pointer(en);
    KernelInfo *mykernel = AcquireKernelInfo(ks);
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
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix"
          " dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    mdr_param = ent_data(en);
    if (mdr_param->nrow != mdr_param->ncol)
    {
      fprintf(stderr, THIS_SOLVER  ": method %s requires 'matrix' real square matrix"
          " dim 3-5 as second argument!\n", method);
      rerror (THIS_SOLVER  ": Cannot continue. Parameters are missing!\n");
    }
    order = mdr_Transpose(mdr_param);

    if (ch != UndefinedChannel)
    {
      MagickConvolveImageChannel(magick_wand[default_wand_idx], ch,
                                 (size_t) MNR(order), MDRPTR(order));
    }
    else
    {
      MagickConvolveImage(magick_wand[default_wand_idx], (size_t) MNR(order), MDRPTR(order));
    }

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
          if (ch != UndefinedChannel)
          {
            func_ch(magick_wand[default_wand_idx], ch);
          }
          else
          {
            func(magick_wand[default_wand_idx]);
          }
        }
        break;

      case 1:
        if (!have_ch)
        {
          func(magick_wand[default_wand_idx], params[0]);
        }
        else
        {
          if (ch != UndefinedChannel)
          {
            func_ch(magick_wand[default_wand_idx], ch, params[0]);
          }
          else
          {
            func(magick_wand[default_wand_idx], params[0]);
          }
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
          if (ch != UndefinedChannel)
          {
            func_ch(magick_wand[default_wand_idx], ch,params[0], params[1]);
          }
          else
          {
            func(magick_wand[default_wand_idx], params[0], params[1]);
          }
        }
        break;

      case 3:
        // three parameter functions
        if (!have_ch)
        {
          func(magick_wand[default_wand_idx], params[0], params[1],params[2]);
        }
        else
        {
          if (ch != UndefinedChannel)
          {
            func_ch(magick_wand[default_wand_idx],ch,params[0],params[1],params[2]);
          }
          else
          {
            func(magick_wand[default_wand_idx], params[0], params[1],params[2]);
          }
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
  ChannelType ch=AllChannels;

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
    ch = string2channel(c);
  }

  // if wand j exists destroy it
  if (magick_wand[j])
  {
    magick_wand[j] = DestroyMagickWand(magick_wand[j]);
    magick_wand[j] = 0;
  }

  magick_wand[j] = MagickCombineImages(magick_wand[default_wand_idx], ch);

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
  ChannelType ch=UndefinedChannel;
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

    metric = string2metric(c);
  }

  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3) != MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": Second argument 'channel' must be integer scalar!");
    c = class_char_pointer(e3);
    ch = string2channel(c);
  }

  if (ch != UndefinedChannel)
  {
    magick_wand[n] = MagickCompareImageChannels(magick_wand[default_wand_idx], magick_wand[j],
        ch, metric, &distortion);
  }
  else
  {
    magick_wand[n] = MagickCompareImages(magick_wand[default_wand_idx], magick_wand[j],
                                         metric, &distortion);
  }

  if (magick_wand[n])
  {
    default_wand_idx = n;

    // add distortion to the artifacts of default image
    sprintf(string_buff, "%g", distortion);
    MagickSetImageArtifact(magick_wand[default_wand_idx], "distortion", string_buff);
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
Ent * ent_im_wand_composite (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *x=0;
  int n;
  char *c=0;
  ChannelType ch=UndefinedChannel;
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
    compose = string2compositeoperator(c);
  }

  if (nargs >= 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": Second argument 'offset' must be two-column vector!");

    x = ent_data(e2);
    if (SIZE(x) != 2)
      rerror(THIS_SOLVER ": Second argument 'offset' must be two-column vector!");
  }

  if (nargs >= 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3) != MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": Third argument 'channel' must be integer scalar!");
    c = class_char_pointer(e3);
    ch = string2channel(c);
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

  if (ch != UndefinedChannel)
  {
    MagickCompositeImageChannel(magick_wand[n], ch, magick_wand[default_wand_idx], compose,
                                (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));
  }
  else
  {
    MagickCompositeImage(magick_wand[n], magick_wand[default_wand_idx], compose,
                         (size_t) MdrV0(x,0), (size_t) MdrV0(x,1));
  }

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
// combine image with an image containing color lookup table
//
#undef THIS_SOLVER
#define THIS_SOLVER RLAB_NAME_EMBED_IM ".clut"
Ent * ent_im_wand_clut (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  int j;
  char *c=0;
  ChannelType ch=UndefinedChannel;

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
    ch = string2channel(c);
  }

  if (ch != UndefinedChannel)
  {
    MagickClutImageChannel(magick_wand[default_wand_idx], ch, magick_wand[j]);
  }
  else
  {
    MagickClutImage(magick_wand[default_wand_idx], magick_wand[j]);
  }

  ent_Clean(e1);
  ent_Clean(e2);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
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
  int is_eval=0, is_func=0, is_fx=0, len=0;

  MagickEvaluateOperator op=UndefinedEvaluateOperator;
  MagickFunction function=UndefinedFunction;
  ChannelType ch=UndefinedChannel;

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
  is_eval = 1;
  op = string2op(expr);

  if (op == UndefinedEvaluateOperator)
  {
    is_eval = 0;
    function = string2func(expr);
    is_func = 1;
  }
  else
  {
    // assume it is expression and just pass it along
    is_func = 0;
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
      ch = string2channel(c);
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
      len = SIZE(params);
    }
  }

  if (is_eval)
  {
    if (ch != UndefinedChannel)
    {
      //
      // call evaluate
      //
      MagickEvaluateImageChannel(magick_wand[default_wand_idx], ch, op, MdrV0(params,0));
    }
    else
    {
      MagickEvaluateImage(magick_wand[default_wand_idx], op, MdrV0(params,0));
    }
  }
  else if (is_fx)
  {
    if (ch != UndefinedChannel)
    {
      //
      // call fx
      //
      MagickFxImageChannel(magick_wand[default_wand_idx], ch, expr);
    }
    else
    {
      MagickFxImage(magick_wand[default_wand_idx], expr);
    }
  }
  else if (is_func)
  {
    if (ch != UndefinedChannel)
    {
      MagickFunctionImageChannel(magick_wand[default_wand_idx], ch,
                                 function, (size_t) len, MDRPTR(params));
    }
    else
    {
      MagickFunctionImage(magick_wand[default_wand_idx], function, (size_t) len, MDRPTR(params));
    }
  }

  ent_Clean(e1);
  ent_Clean(e2);

  rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (default_wand_idx + 1);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

