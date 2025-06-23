# This document describes how to install rlabplus on raspbian/debian linux systems.
# Items listed in [..] might be part of default installation (e.g., noobs).
# Following packages are needed:

apt-get install \
  gawk  gfortran  libcurl4-openssl-dev \
  [libblas-common] libblas-dev liblapack-dev \
  libgsl-dev libhdf5-dev libhdf5-serial-dev \
  libreadline-dev libncurses-dev libgphoto2-dev \
  imagemagick libmagick-dev libmagickcore-6-headers libmagickwand-dev \
  libglpk-dev libarpack2-dev libglib2.0-dev libwebp-dev \
  libx11-dev libxt-dev libsqlite3-dev libssl-dev  \
  bison flex

# these plotting programs may be useful when working with rlab
apt-get install \
  xfig gnuplot xmgr

# then do just for a good measure
make clean

# start with
./rconfigure3_debian_rpi

# then decide on the parser
make scanner 

# then continue to build parts, in no particular order,
make gc clibs flibs

# this one goes last
make rlab

# install/uninstall it
sudo make (un)install

