#-------------------------------------------------------------------------------
# UMFPACK makefile (for GNU make or original make)
#-------------------------------------------------------------------------------

# UMFPACK requires the AMD package to be in ../AMD

# compile everything (except hb, fortran, and fortran64), including AMD and the
# MATLAB mexFunctions
all:
	( cd ../AMD ; make )
	( cd Source ; make )

# compile just the C-callable libraries and demo programs (not mexFunctions)
lib:
	( cd ../AMD ; make lib )
	( cd Source ; make )

# remove object files, but keep the compiled programs and library archives
cleanlib:
	( cd ../AMD ; make clean )
	( cd Source ; make clean )
	( cd Doc    ; make clean )

# clean, and then remove compiled programs and library archives
clean:
	( cd ../AMD ; make purge )
	( cd Source ; make purge )
	( cd Doc    ; make purge )

# create PDF documents for the original distribution
doc:
	( cd ../AMD ; make doc )
	( cd Doc    ; make )

# get ready for distribution
dist: purge
	( cd ../AMD ; make dist )
	( cd Doc    ; make )
