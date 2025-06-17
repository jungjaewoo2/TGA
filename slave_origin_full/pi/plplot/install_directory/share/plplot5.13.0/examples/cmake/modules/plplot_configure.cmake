# -*- mode: cmake -*-
# This file used to configure CMake variables that are used to build the
# installed examples.
# Identify that this is the build system for the installed examples.
set(BUILD_TEST ON)
set(CORE_BUILD OFF)

# =======================================================================
# Useful functions.
# =======================================================================
include(plplot_functions)

# =======================================================================
# pkg-config support as well as macros to put link flags in standard
# *.pc (pkg-config) form as well as standard fullpath form used by cmake.
# PKG_CONFIG_EXECUTABLE can be used to discover whether pkg-config was
# found or not.
# =======================================================================
include(pkg-config)

set(INCLUDE_DIR "/home/pi/plplot/install_directory/include/plplot")
set(CMAKE_PREFIX_PATH "/home/pi/plplot/install_directory")
set(MATH_LIB /usr/lib/arm-linux-gnueabihf/libm.so)
set(SH_EXECUTABLE /bin/bash)

set(BUILD_SHARED_LIBS ON)

set(language_info_LIST ada:ada:a;c:c:c;cxx:c++:cxx;d:d:d;fortran:fortran:f;java:java:j;lua:lua:lua;ocaml:ocaml:ocaml;octave:octave:o;pdl:perl:pdl;python:python:p;tcl:tcl:t;plrender:c:plm)

# This list of set commands must be consistent with language_info_LIST
# set in examples/CMakeLists.txt and configured just above.  Thus,
# some consistent maintenance of that list and the current list of set
# commands must be done on the rare occasions when a language is added
# to PLplot.  Cannot make the list of set commands automatically with
# a foreach loop because explicit configuration of each name is required.
set(ENABLE_ada OFF)
set(ENABLE_c ON)
set(ENABLE_cxx ON)
set(ENABLE_d OFF)
set(ENABLE_fortran ON)
set(ENABLE_java OFF)
set(ENABLE_lua OFF)
set(ENABLE_ocaml OFF)
set(ENABLE_octave OFF)
set(ENABLE_pdl OFF)
set(ENABLE_python OFF)
set(ENABLE_pyqt4 OFF)
set(ENABLE_pyqt5 OFF)
set(ENABLE_tcl OFF)

# Interactive bindings
set(ENABLE_tk OFF)
set(ENABLE_tkX OFF)
set(OLD_WXWIDGETS )
set(wxdemo_name )
set(ENABLE_wxwidgets OFF)
set(ENABLE_qt OFF)

# Check for all officially supported CMake compilers that we use with
# soft-landing support.
include(language_support)

# First field is lower-case language string, second field is
# mixed-case language string used for official CMake variables, third
# field is language string to be used in output messages, and final
# field is whether it is a fatal error if this language is missing/not
# working.
set(COMPILED_LANGUAGE_INFO_LIST
"c:C:C:ON"
"ada:Ada:Ada:OFF"
"cxx:CXX:C++:OFF"
"d:D:D:OFF"
"fortran:Fortran:Fortran:OFF"
"java:Java:Java:OFF"
)

# C++ is a special language case that needs to be enabled if at all
# possible regardless of ENABLE_CXX (see cmake/modules/c++.cmake.)
# Hence the (language STREQUAL "cxx" OR ENABLE_${language}) logic below.

foreach(COMPILED_LANGUAGE_INFO ${COMPILED_LANGUAGE_INFO_LIST})
  string(REGEX REPLACE "^(.*):.*:.*:.*$" "\\1" language ${COMPILED_LANGUAGE_INFO})
  string(REGEX REPLACE "^.*:(.*):.*:.*$" "\\1" Language ${COMPILED_LANGUAGE_INFO})
  string(REGEX REPLACE "^.*:.*:(.*):.*$" "\\1" language_output ${COMPILED_LANGUAGE_INFO})
  string(REGEX REPLACE "^.*:.*:.*:(.*)$" "\\1" fatal_language ${COMPILED_LANGUAGE_INFO})
  if((language STREQUAL "cxx" OR ENABLE_${language}) AND NOT PLPLOT_${Language}_COMPILER_WORKS)
    workaround_9220(${Language} PLPLOT_${Language}_COMPILER_WORKS)
    if(NOT PLPLOT_${Language}_COMPILER_WORKS)
      if(fatal_language)
	message(FATAL_ERROR "${language_output} compiler absolutely required to build this project.")
      endif(fatal_language)
      message(STATUS "WARNING: no working ${language_output} compiler so disabling ${language} examples.")
      set(ENABLE_${language} OFF CACHE BOOL "Enable ${language} bindings" FORCE)
    endif(NOT PLPLOT_${Language}_COMPILER_WORKS)
  endif((language STREQUAL "cxx" OR ENABLE_${language}) AND NOT PLPLOT_${Language}_COMPILER_WORKS)

  if(language STREQUAL "cxx" OR ENABLE_${language})
    # Find and check ${Language} compiler
    enable_language(${Language} OPTIONAL)
    if(NOT CMAKE_${Language}_COMPILER_WORKS)
      if(fatal_language)
	message(FATAL_ERROR "${language_output} compiler absolutely required to build this project.")
      endif(fatal_language)
      message(STATUS "WARNING: no working ${language_output} compiler so disabling ${language} examples.")
      set(ENABLE_${language} OFF CACHE BOOL "Enable ${language} bindings" FORCE)
    endif(NOT CMAKE_${Language}_COMPILER_WORKS)
  endif(language STREQUAL "cxx" OR ENABLE_${language})
endforeach(COMPILED_LANGUAGE_INFO ${COMPILED_LANGUAGE_INFO_LIST})

if(ENABLE_ada)
  set(ADA_INCLUDE_DIR "")
  set(ADA_LIB_DIR "")
endif(ENABLE_ada)

if(ENABLE_fortran)
  set(FORTRAN_MOD_DIR "/home/pi/plplot/install_directory/lib/fortran/modules/plplot")
endif(ENABLE_fortran)

if(ENABLE_java)
  find_package(JNI)
  set(JAR_DIR "")
  set(LIB_DIR "/home/pi/plplot/install_directory/lib")
endif(ENABLE_java)

if(ENABLE_ocaml)
  set(OCAMLC OCAMLC-NOTFOUND)
  set(OCAMLOPT )
  set(OCAMLFIND )
  set(OCAML_HAS_CAIRO OFF)
  set(OCAML_HAS_GTK OFF)
  set(CMAKE_INSTALL_LIBDIR "/home/pi/plplot/install_directory/lib")
  set(OCAML_INSTALL_DIR "")
endif(ENABLE_ocaml)

if(ENABLE_tcl OR ENABLE_tk)
  set(TCL_TCLSH "TCL_TCLSH-NOTFOUND")
  set(MKTCLINDEX ${CMAKE_SOURCE_DIR}/tcl/mktclIndex)
  set(MKTCLINDEX_ARGS )
  # Needed to configure tclsh_standard_examples and wish_standard_examples
  set(DATA_DIR "/home/pi/plplot/install_directory/share/plplot5.13.0")
endif(ENABLE_tcl OR ENABLE_tk)

if(ENABLE_tk)
  set(ENABLE_itk OFF)
  set(ENABLE_itkX OFF)
  set(TCL_INCLUDE_PATH "TCL_INCLUDE_PATH-NOTFOUND")
  set(TK_INCLUDE_PATH "TK_INCLUDE_PATH-NOTFOUND")
  set(TCL_LIBRARY "/usr/lib/arm-linux-gnueabihf/libtcl8.6.so")
  set(TK_LIBRARY "/usr/lib/arm-linux-gnueabihf/libtk8.6.so")
endif(ENABLE_tk)

if(ENABLE_wxwidgets)
  # Use identical compile and link flags to build wxwidgets application.
  set(wxwidgets_COMPILE_FLAGS "")
  set(wxwidgets_LINK_FLAGS )
endif(ENABLE_wxwidgets)

if(ENABLE_qt)
  set(ANY_QT_DEVICE OFF)
  # Qt5 support?
  set(PLPLOT_USE_QT5 ON)
  if(PLPLOT_USE_QT5)
    set(CORE_Qt5_VERSION_MAJOR 5)
    set(CORE_Qt5_VERSION_MINOR 7)
    set(CORE_Qt5_VERSION_PATCH 1)
  else(PLPLOT_USE_QT5)
    set(CORE_QT_VERSION_MAJOR )
    set(CORE_QT_VERSION_MINOR )
    set(CORE_QT_VERSION_PATCH )
  endif(PLPLOT_USE_QT5)
  set(PLPLOT_AUTOMOC_MOC_OPTIONS )
endif(ENABLE_qt)

set(ENABLE_DYNDRIVERS ON)

set(FILE_DEVICES_LIST ps:ps:OFF;svg:svg:ON;xfig:xfig:ON)
# This list of set commands must be consistent with the above list which
# implies a maintenance issue.
set(PLD_pdfcairo OFF)
set(PLD_pngcairo OFF)
set(PLD_pscairo OFF)
set(PLD_epscairo OFF)
set(PLD_svgcairo OFF)
set(PLD_cgm OFF)
set(PLD_epsqt OFF)
set(PLD_pdfqt OFF)
set(PLD_bmpqt OFF)
set(PLD_jpgqt OFF)
set(PLD_pngqt OFF)
set(PLD_ppmqt OFF)
set(PLD_tiffqt OFF)
set(PLD_svgqt OFF)
set(PLD_gif OFF)
set(PLD_jpeg OFF)
set(PLD_png OFF)
set(PLD_pdf OFF)
set(PLD_ps ON)
set(PLD_pstex OFF)
set(PLD_psttf OFF)
set(PLD_svg ON)
set(PLD_wxpng OFF)
set(PLD_xfig ON)

# These only used for testing, but not core build
set(PLD_psc ON)
set(PLD_psttfc OFF)

set(INTERACTIVE_DEVICES_LIST xwin:xwin)
# Set interactive devices that are needed to build certain examples.  This
# list should be consistent with INTERACTIVE_DEVICES_LIST which implies a
# maintenance issue.

set(PLD_aqt OFF)
set(PLD_ntk OFF)
set(PLD_qtwidget OFF)
set(PLD_tk OFF)
set(PLD_wincairo OFF)
set(PLD_wingcc OFF)
set(PLD_wxwidgets OFF)
set(PLD_xcairo OFF)
set(PLD_xwin ON)
set(PLD_wingdi OFF)

# Set external devices that are needed to build certain examples.
set(PLD_extcairo OFF)
set(PLD_extqt OFF)

# Needed to check if diff and tail on system for test_diff_psc target.

set(DIFF_EXECUTABLE /usr/bin/diff)
set(TAIL_EXECUTABLE /usr/bin/tail)
set(CMP_EXECUTABLE /usr/bin/cmp)

# ocaml examples build needs this.
set(USE_RPATH ON)

# WIN32 and CYGWIN derived variables
set(WIN32_AND_NOT_CYGWIN )
set(WIN32_OR_CYGWIN OFF)

# Needed for Tcl/Tk linking decisions.
set(tcltk_in_plplot_library OFF)

# Needed for running Python scripts if PYTHON_EXECUTABLE is different
# from /usr/bin/env python (which often happens if both python2 and python3
# are installed).
set(PYTHON_EXECUTABLE )
# Needed for Python testing
set(PYTHON_VERSION )
