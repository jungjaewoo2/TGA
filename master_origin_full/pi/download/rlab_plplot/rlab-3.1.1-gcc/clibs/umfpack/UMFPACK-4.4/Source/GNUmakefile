#-------------------------------------------------------------------------------
# UMFPACK Makefile for compiling on Unix systems (for GNU Make)
#-------------------------------------------------------------------------------

all: ../Lib/libumfpack.a

include ../Make/Make.include

C = $(CC) $(CFLAGS) $(CONFIG) -I../Include -I../../AMD/Include -I../../AMD/Source

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

# non-user-callable umf_*.[ch] files:
UMFCH = umf_assemble umf_blas3_update umf_build_tuples umf_create_element \
	umf_dump umf_extend_front umf_garbage_collection umf_get_memory \
	umf_init_front umf_kernel umf_kernel_init umf_kernel_wrapup \
	umf_local_search umf_lsolve umf_ltsolve umf_mem_alloc_element \
	umf_mem_alloc_head_block umf_mem_alloc_tail_block \
	umf_mem_free_tail_block umf_mem_init_memoryspace \
	umf_report_vector umf_row_search umf_scale_column \
	umf_set_stats umf_solve umf_symbolic_usage umf_transpose \
	umf_tuple_lengths umf_usolve umf_utsolve umf_valid_numeric \
	umf_valid_symbolic umf_grow_front umf_start_front umf_2by2 \
	umf_store_lu umf_scale

# non-user-callable umf_*.[ch] files, int/long versions only (no real/complex):
UMFINT = umf_analyze umf_apply_order umf_colamd umf_free umf_fsize \
	umf_is_permutation umf_malloc umf_realloc umf_report_perm \
	umf_singletons

# non-user-callable, created from umf_ltsolve.c, umf_utsolve.c,
# umf_triplet.c, and umf_assemble.c , with int/long and real/complex versions:
UMF_CREATED = umf_lhsolve umf_uhsolve umf_triplet_map_nox \
	umf_triplet_nomap_x umf_triplet_nomap_nox umf_triplet_map_x \
	umf_assemble_fixq umf_store_lu_drop

# non-user-callable, int/long and real/complex versions:
UMF = $(UMF_CREATED) $(UMFCH)

# user-callable umfpack_*.[ch] files (int/long and real/complex):
UMFPACK = umfpack_col_to_triplet umfpack_defaults umfpack_free_numeric \
	umfpack_free_symbolic umfpack_get_numeric umfpack_get_lunz \
	umfpack_get_symbolic umfpack_get_determinant umfpack_numeric \
	umfpack_qsymbolic umfpack_report_control umfpack_report_info \
	umfpack_report_matrix umfpack_report_numeric umfpack_report_perm \
	umfpack_report_status umfpack_report_symbolic umfpack_report_triplet \
	umfpack_report_vector umfpack_solve umfpack_symbolic \
	umfpack_transpose umfpack_triplet_to_col umfpack_scale \
	umfpack_load_numeric umfpack_save_numeric \
	umfpack_load_symbolic umfpack_save_symbolic

# user-callable, created from umfpack_solve.c (umfpack_wsolve.h exists, though):
# with int/long and real/complex versions:
UMFPACKW = umfpack_wsolve

USER = $(UMFPACKW) $(UMFPACK)

# user-callable, only one version for int/long, real/complex, *.[ch] files:
GENERIC = umfpack_timer umfpack_tictoc

#-------------------------------------------------------------------------------
# include files:
#-------------------------------------------------------------------------------

INC = ../Include/umfpack.h \
	umf_config.h umf_version.h umf_internal.h umf_triplet.h \
	$(addsuffix .h,$(UMFCH)) \
	$(addsuffix .h,$(UMFINT)) \
	$(addprefix ../Include/, $(addsuffix .h,$(USER))) \
	$(addprefix ../Include/, $(addsuffix .h,$(GENERIC))) \
	../../AMD/Source/amd_internal.h ../../AMD/Include/amd.h

#-------------------------------------------------------------------------------
# object files for each version
#-------------------------------------------------------------------------------

DI = $(addsuffix .o, $(subst umf_,umf_di_,$(UMF)) $(subst umfpack_,umfpack_di_,$(USER)))
DL = $(addsuffix .o, $(subst umf_,umf_dl_,$(UMF)) $(subst umfpack_,umfpack_dl_,$(USER)))
ZI = $(addsuffix .o, $(subst umf_,umf_zi_,$(UMF)) $(subst umfpack_,umfpack_zi_,$(USER)))
ZL = $(addsuffix .o, $(subst umf_,umf_zl_,$(UMF)) $(subst umfpack_,umfpack_zl_,$(USER)))
I  = $(addsuffix .o, $(subst umf_,umf_i_,$(UMFINT)))
L  = $(addsuffix .o, $(subst umf_,umf_l_,$(UMFINT)))
GN = $(addsuffix .o, $(subst umfpack_,umfpack_gn_,$(GENERIC)))

#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

umf_i_%.o: umf_%.c $(INC)
	$(C) -DDINT -c $< -o $@

umf_l_%.o: umf_%.c $(INC)
	$(C) -DDLONG -c $< -o $@

#-------------------------------------------------------------------------------
# compile each routine in the DI version
#-------------------------------------------------------------------------------

umf_di_%.o: umf_%.c $(INC)
	$(C) -DDINT -c $< -o $@

umf_di_%hsolve.o: umf_%tsolve.c $(INC)
	$(C) -DDINT -DCONJUGATE_SOLVE -c $< -o $@

umf_di_triplet_map_x.o: umf_triplet.c $(INC)
	$(C) -DDINT -DDO_MAP -DDO_VALUES -c $< -o $@

umf_di_triplet_map_nox.o: umf_triplet.c $(INC)
	$(C) -DDINT -DDO_MAP -c $< -o $@

umf_di_triplet_nomap_x.o: umf_triplet.c $(INC)
	$(C) -DDINT -DDO_VALUES -c $< -o $@

umf_di_triplet_nomap_nox.o: umf_triplet.c $(INC)
	$(C) -DDINT -c $< -o $@

umf_di_assemble_fixq.o: umf_assemble.c $(INC)
	$(C) -DDINT -DFIXQ -c $< -o $@

umf_di_store_lu_drop.o: umf_store_lu.c $(INC)
	$(C) -DDINT -DDROP -c $< -o $@

umfpack_di_wsolve.o: umfpack_solve.c $(INC)
	$(C) -DDINT -DWSOLVE -c $< -o $@

umfpack_di_%.o: umfpack_%.c $(INC)
	$(C) -DDINT -c $< -o $@

#-------------------------------------------------------------------------------
# compile each routine in the DL version
#-------------------------------------------------------------------------------

umf_dl_%.o: umf_%.c $(INC)
	$(C) -DDLONG -c $< -o $@

umf_dl_%hsolve.o: umf_%tsolve.c $(INC)
	$(C) -DDLONG -DCONJUGATE_SOLVE -c $< -o $@

umf_dl_triplet_map_x.o: umf_triplet.c $(INC)
	$(C) -DDLONG -DDO_MAP -DDO_VALUES -c $< -o $@

umf_dl_triplet_map_nox.o: umf_triplet.c $(INC)
	$(C) -DDLONG -DDO_MAP -c $< -o $@

umf_dl_triplet_nomap_x.o: umf_triplet.c $(INC)
	$(C) -DDLONG -DDO_VALUES -c $< -o $@

umf_dl_triplet_nomap_nox.o: umf_triplet.c $(INC)
	$(C) -DDLONG -c $< -o $@

umf_dl_assemble_fixq.o: umf_assemble.c $(INC)
	$(C) -DDLONG -DFIXQ -c $< -o $@

umf_dl_store_lu_drop.o: umf_store_lu.c $(INC)
	$(C) -DDLONG -DDROP -c $< -o $@

umfpack_dl_wsolve.o: umfpack_solve.c $(INC)
	$(C) -DDLONG -DWSOLVE -c $< -o $@

umfpack_dl_%.o: umfpack_%.c $(INC)
	$(C) -DDLONG -c $< -o $@

#-------------------------------------------------------------------------------
# compile each routine in the ZI version
#-------------------------------------------------------------------------------

umf_zi_%.o: umf_%.c $(INC)
	$(C) -DZINT -c $< -o $@

umf_zi_%hsolve.o: umf_%tsolve.c $(INC)
	$(C) -DZINT -DCONJUGATE_SOLVE -c $< -o $@

umf_zi_triplet_map_x.o: umf_triplet.c $(INC)
	$(C) -DZINT -DDO_MAP -DDO_VALUES -c $< -o $@

umf_zi_triplet_map_nox.o: umf_triplet.c $(INC)
	$(C) -DZINT -DDO_MAP -c $< -o $@

umf_zi_triplet_nomap_x.o: umf_triplet.c $(INC)
	$(C) -DZINT -DDO_VALUES -c $< -o $@

umf_zi_triplet_nomap_nox.o: umf_triplet.c $(INC)
	$(C) -DZINT -c $< -o $@

umf_zi_assemble_fixq.o: umf_assemble.c $(INC)
	$(C) -DZINT -DFIXQ -c $< -o $@

umf_zi_store_lu_drop.o: umf_store_lu.c $(INC)
	$(C) -DZINT -DDROP -c $< -o $@

umfpack_zi_wsolve.o: umfpack_solve.c $(INC)
	$(C) -DZINT -DWSOLVE -c $< -o $@

umfpack_zi_%.o: umfpack_%.c $(INC)
	$(C) -DZINT -c $< -o $@

#-------------------------------------------------------------------------------
# compile each routine in the ZL version
#-------------------------------------------------------------------------------

umf_zl_%.o: umf_%.c $(INC)
	$(C) -DZLONG -c $< -o $@

umf_zl_%hsolve.o: umf_%tsolve.c $(INC)
	$(C) -DZLONG -DCONJUGATE_SOLVE -c $< -o $@

umf_zl_triplet_map_x.o: umf_triplet.c $(INC)
	$(C) -DZLONG -DDO_MAP -DDO_VALUES -c $< -o $@

umf_zl_triplet_map_nox.o: umf_triplet.c $(INC)
	$(C) -DZLONG -DDO_MAP -c $< -o $@

umf_zl_triplet_nomap_x.o: umf_triplet.c $(INC)
	$(C) -DZLONG -DDO_VALUES -c $< -o $@

umf_zl_triplet_nomap_nox.o: umf_triplet.c $(INC)
	$(C) -DZLONG -c $< -o $@

umf_zl_assemble_fixq.o: umf_assemble.c $(INC)
	$(C) -DZLONG -DFIXQ -c $< -o $@

umf_zl_store_lu_drop.o: umf_store_lu.c $(INC)
	$(C) -DZLONG -DDROP -c $< -o $@

umfpack_zl_wsolve.o: umfpack_solve.c $(INC)
	$(C) -DZLONG -DWSOLVE -c $< -o $@

umfpack_zl_%.o: umfpack_%.c $(INC)
	$(C) -DZLONG -c $< -o $@

#-------------------------------------------------------------------------------
# Create the generic routines (GN) using a generic rule
#-------------------------------------------------------------------------------

umfpack_gn_%.o: umfpack_%.c $(INC)
	$(C) -c $< -o $@

#-------------------------------------------------------------------------------
# Create the libumfpack.a library
#-------------------------------------------------------------------------------

../Lib/libumfpack.a: $(I) $(L) $(GN) $(DI) $(DL) $(ZI) $(ZL)
	$(AR) cr ../Lib/libumfpack.a $^
	- $(RANLIB) ../Lib/libumfpack.a

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

purge: clean
	- $(RM) ../Lib/libumfpack.a
