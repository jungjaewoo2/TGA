#
# Simple Makefile for Rlab documentation.
# The doc sources are all in SGML (linuxdoc-sgml)
# and are translated to various formats...
#

SGMLFILES =	ref.sgml intro.sgml data-classes.sgml io.sgml \
		control-flow.sgml ufunctions.sgml pref.sgml tutorial.sgml

HELPSRC = 	help-source/abs.sgml	help-source/acos.sgml \
		help-source/all.sgml	help-source/any.sgml \
		help-source/ascii.sh	help-source/asin.sgml \
		help-source/atan.sgml	help-source/atan2.sgml \
		help-source/backsub.sgml help-source/balance.sgml \
		help-source/ceil.sgml	help-source/chol.sgml \
		help-source/class.sgml	help-source/clear.sgml \
		help-source/close.sgml	help-source/conj.sgml \
		help-source/cos.sgml	help-source/cumprod.sgml \
		help-source/cumsum.sgml	help-source/det.sgml \
		help-source/diag.sgml	help-source/dlopen.sgml \
		help-source/eig.sgml	help-source/entinfo.sgml \
		help-source/error.sgml	help-source/eval.sgml \
		help-source/exist.sgml	help-source/exp.sgml \
		help-source/factor.sgml	help-source/fft.sgml \
		help-source/filter.sgml	help-source/find.sgml \
		help-source/finite.sgml	help-source/floor.sgml \
		help-source/format.sgml	help-source/fprintf.sgml \
		help-source/fread.sgml	help-source/frexp.sgml \
		help-source/fseek.sgml	help-source/full.sgml \
		help-source/getenv.sgml	help-source/getline.sgml \
		help-source/help.sgml	help-source/hess.sgml \
		help-source/ifft.sgml	help-source/imag.sgml \
		help-source/inf.sgml	help-source/int.sgml \
		help-source/isinf.sgml	help-source/isnan.sgml \
		help-source/issymm.sgml	help-source/ldexp.sgml \
		help-source/length.sgml \
		help-source/load.sgml	help-source/log.sgml \
		help-source/log10.sgml	help-source/logb.sgml \
		help-source/max.sgml	help-source/maxi.sgml \
		help-source/members.sgml	help-source/min.sgml \
		help-source/mini.sgml	help-source/mnorm.sgml \
		help-source/mod.sgml	help-source/nan.sgml \
		help-source/ode.sgml	help-source/ones.sgml \
		help-source/open.sgml	help-source/printf.sgml \
		help-source/prod.sgml	help-source/putenv.sgml \
		help-source/qr.sgml	help-source/quit.sgml \
		help-source/rand.sgml	help-source/rcond.sgml \
		help-source/read.sgml	help-source/read_ascii.sgml \
		help-source/readm.sgml \
		help-source/real.sgml	help-source/require.sgml \
		help-source/reshape.sgml	help-source/rfile.sgml \
		help-source/round.sgml	help-source/schur.sgml \
		help-source/sgml2help.sh	help-source/sgml2page.sh \
		help-source/sign.sgml	help-source/sin.sgml \
		help-source/size.sgml	help-source/sizeof.sgml \
		help-source/solve.sgml	help-source/sort.sgml \
		help-source/sparse.sgml	help-source/spconvert.sgml \
		help-source/sprintf.sgml	help-source/sqrt.sgml \
		help-source/srand.sgml	help-source/strsplt.sgml \
		help-source/strtod.sgml	help-source/sum.sgml \
		help-source/svd.sgml	help-source/system.sgml \
		help-source/tan.sgml	help-source/tic.sgml \
		help-source/tmpnam.sgml	help-source/toc.sgml \
		help-source/type.sgml	help-source/write.sgml \
		help-source/write_ascii.sgml \
		help-source/writem.sgml	help-source/zeros.sgml \
		help-source/vpnorm.sgml help-source/spwrite.sgml \
		help-source/fwrite.sgml	help-source/strtol.sgml

#
# Run include.awk twice. The second time to pick up any includes
# in the first set of include files.
#

all:		rlab-ref.sgml rlab-ref.dvi htmlfiles helpfiles

rlab-ref.sgml:	$(SGMLFILES) $(HELPSRC)
		mawk -f include.awk < ref.sgml > rlab-ref.sgml
		mawk -f include.awk < rlab-ref.sgml > tmp.sgml
		mv tmp.sgml rlab-ref.sgml

rlab-ref.dvi:	rlab-ref.sgml	
		sgml2latex -d rlab-ref.sgml

htmlfiles:	rlab-ref.sgml $(SGMLFILES)
		sgml2html rlab-ref.sgml
		mv *.html html/

.PHONY:		htmlfiles

helpfiles:	
		./make-help-files.sh

.PHONY:		helpfiles

