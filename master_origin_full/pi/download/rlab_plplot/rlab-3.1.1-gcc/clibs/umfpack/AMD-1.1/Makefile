#------------------------------------------------------------------------------
# AMD Makefile (for GNU Make or original make)
#------------------------------------------------------------------------------

# Compile everything, including the C-callable routine.
all:
	( cd Source ; make )

# compile just the C-callable libraries and demo programs (not mexFunctions)
lib:
	( cd Source ; make )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd Source ; make clean )
	( cd Doc    ; make clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd Source ; make purge )
	( cd Doc    ; make purge )

# create PDF documents for the original distribution
doc:
	( cd Doc    ; make )

# get ready for distribution
dist: purge
	( cd Doc    ; make )
