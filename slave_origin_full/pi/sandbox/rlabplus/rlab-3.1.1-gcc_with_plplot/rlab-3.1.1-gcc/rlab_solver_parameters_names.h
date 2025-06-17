// Copyright (C) 2003-2013 Marijan Kostrun
//   part of rlabplus project on rlabplus.sourceforge.net
//
// solver parameters names
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// See the file ./COPYING
// **********************************************************************

#ifndef RLAB_PARAMETER_NAMES_ETC
#define RLAB_PARAMETER_NAMES_ETC


// The purpose of this file is to provide the names for the
// solver parameters in rlabplus. This is because rlab has not
// reached yet convergence regarding what should be certain parameters
// called, on one hand side, and on the other, to provide an overview
// for the developer on what names have been chosen so far. The latter
// is in attempt to bring some uniformity to naming procedures and practices.

// what algorithm to use for sorting
#define RLAB_SORT_QUICK 0
#define RLAB_SORT_HEAP 1

// how to sort
#define RLAB_SORT_NAN_INPLACE   0
#define RLAB_SORT_NAN_ONBOTTOM -1
#define RLAB_SORT_NAN_ONTOP     1

// Add this line to the 'define/include' section of the c-file
// naming convention for the solver parameters
#define RLAB_STATUS_SUCCESS 0
#define RLAB_STATUS_FAILURE 1
#define RLAB_STATUS_ERROR   (-1)

#define RLAB_SPRINTF_NAN "NaN"

#define RLAB_SYSVAR_CONFIG          "_rlab_config"
#define RLAB_SYSVAR_CONFIG_VERSION  "ver"
#define RLAB_SYSVAR_CONFIG_PLOT     "plot_support"

#define RLAB_MEMBER_NROW            "nr"
#define RLAB_MEMBER_NCOL            "nc"
#define RLAB_MEMBER_SIZE            "n"
#define RLAB_MEMBER_SIZE_NON_ZERO   "nnz"
#define RLAB_MEMBER_CLASS           "class"
#define RLAB_MEMBER_CLASS_NUM       "num"
#define RLAB_MEMBER_CLASS_F         "function"
#define RLAB_MEMBER_CLASS_LIST      "list"
#define RLAB_MEMBER_TYPE            "type"
#define RLAB_MEMBER_TYPE_REAL       "real"
#define RLAB_MEMBER_TYPE_COMPLEX    "complex"
#define RLAB_MEMBER_TYPE_INT32      "int"
#define RLAB_MEMBER_TYPE_F_BUILTIN  "builtin"
#define RLAB_MEMBER_TYPE_F_USER     "user"
#define RLAB_MEMBER_STORAGE "storage"
#define RLAB_STORAGE_DENSE  "dense"
#define RLAB_STORAGE_SPARSE "sparse"

#define RLAB_CLASS_STRING   "string"
#define RLAB_CLASS_MDE      "cell"
#define RLAB_TYPE_MDE       "cell"

#define RLAB_ENTINFO_ADDR     "addr"
#define RLAB_ENTINFO_DATAPTR  "data_ptr"
#define RLAB_ENTINFO_REFC     "refc"
#define RLAB_ENTINFO_STAT     "_static"



// generic warning messages for users
#define RLAB_ERROR_ARG1_ONLY                     "One argument allowed !"
#define RLAB_ERROR_ARG1_FUNC_VAR                 "1st argument must be function-variable !"
#define RLAB_ERROR_ARG1_NUM                      "1st argument must be numeric matrix !"
#define RLAB_ERROR_ARG1_VAR                      "1st argument must be variable !"
#define RLAB_ERROR_ARG1_MDR                      "1st argument must be real !"
#define RLAB_ERROR_ARG1_MDR_SCALAR               "1st argument must be real scalar !"
#define RLAB_ERROR_ARG1_MDR_REALPOSITIVESCALAR   "1st argument must be real positive scalar !"
#define RLAB_ERROR_ARG1_MDR_VECTOR               "1st argument must be real vector !"
#define RLAB_ERROR_ARG1_MDR_MATRIX               "1st argument must be real matrix !"
#define RLAB_ERROR_ARG1_MDR_MATRIX_2COL          "1st argument must be two-column real matrix !"
#define RLAB_ERROR_ARG1_MDR_VECTOR_4COL          "1st argument must be four-column real vector !"
#define RLAB_ERROR_ARG1_MDR_INTEGER              "1st argument must be integer !"
#define RLAB_ERROR_ARG1_MDR_INTEGER_MATRIX       "1st argument must be integer matrix !"
#define RLAB_ERROR_ARG1_MDR_INTEGER_VECTOR       "1st argument must be integer vector !"
#define RLAB_ERROR_ARG1_MDC                      "1st argument must be complex !"
#define RLAB_ERROR_ARG1_MDS                      "1st argument must be MDS !"
#define RLAB_ERROR_ARG1_MDS_SCALAR               "1st argument must be string !"
#define RLAB_ERROR_ARG1_MDS_VECTOR               "1st argument must be string vector !"
#define RLAB_ERROR_ARG1_MDS_MATRIX               "1st argument must be string matrix !"
#define RLAB_ERROR_ARG1_MD                       "1st argument must be string or numeric matrix !"
#define RLAB_ERROR_ARG1_TABFUNC                  "1st argument must be two-column real matrix !"
#define RLAB_ERROR_ARG1_MDR_TWOROWS              "1st argument real-matrix must have at least two rows !"
#define RLAB_ERROR_ARG1_HIST                     "1st argument must be histogram/list <<bin;range>> !"
#define RLAB_ERROR_ARG1_HIST_OR_MDR_VECTOR       "1st argument must be real vector or histogram/list <<bin;range>> !"

#define RLAB_ERROR_ARG2_FUNC_VAR                 "2nd argument must be function-variable !"
#define RLAB_ERROR_ARG2_MDR                      "2nd argument must be real !"
#define RLAB_ERROR_ARG2_MDR_INTEGER              "2nd argument must be integer !"
#define RLAB_ERROR_ARG2_MDR_INTEGER_VECTOR       "2nd argument must be integer vector!"
#define RLAB_ERROR_ARG2_MDR_INTEGER_MATRIX       "2nd argument must be integer matrix !"
#define RLAB_ERROR_ARG2_MDR_SCALAR               "2nd argument must be real scalar !"
#define RLAB_ERROR_ARG2_VECTOR                   "2nd argument must be vector !"
#define RLAB_ERROR_ARG2_MDR_VECTOR               "2nd argument must be real vector !"
#define RLAB_ERROR_ARG2_MDR_MATRIX               "2nd argument must be real matrix !"
#define RLAB_ERROR_ARG2_MDS                      "2nd argument must be MDS !"
#define RLAB_ERROR_ARG2_MDS_SCALAR               "2nd argument must be string !"
#define RLAB_ERROR_ARG2_MDS_VECTOR               "2nd argument must be string vector !"
#define RLAB_ERROR_ARG2_MDS_MATRIX               "2nd argument must be string matrix !"
#define RLAB_ERROR_ARG2_STRING_OR_REAL_VECTOR    "2nd argument must be string or real vector !"
#define RLAB_ERROR_ARG2_BTREE                    "2nd argument must be list!"

#define RLAB_ERROR_ARG3_FUNC_VAR                 "3rd argument must be function-variable !"
#define RLAB_ERROR_ARG3_MDR                      "3rd argument must be real !"
#define RLAB_ERROR_ARG3_MDR_SCALAR               "3rd argument must be real scalar !"
#define RLAB_ERROR_ARG3_MDR_POSITIVE_SCALAR      "3rd argument must be scalar greater than 0 !"
#define RLAB_ERROR_ARG3_MDR_INTEGER              "3rd argument must be integer !"
#define RLAB_ERROR_ARG3_MDR_INTEGER_MATRIX       "3rd argument must be integer matrix !"
#define RLAB_ERROR_ARG3_MDR_VECTOR               "3rd argument must be real vector !"
#define RLAB_ERROR_ARG3_MDR_MATRIX               "3rd argument must be real matrix !"
#define RLAB_ERROR_ARG3_MDR_2COLMATRIX           "3rd argument must be two-columns real matrix !"
#define RLAB_ERROR_ARG3_MDS                      "3rd argument must be MDS !"
#define RLAB_ERROR_ARG3_MDS_SCALAR               "3rd argument must be string !"
#define RLAB_ERROR_ARG3_MDS_VECTOR               "3rd argument must be string vector !"
#define RLAB_ERROR_ARG3_MDS_MATRIX               "3rd argument must be string matrix !"
#define RLAB_ERROR_ARG3_BTREE                    "3rd argument must be list!"

#define RLAB_ERROR_ARG4_FUNC_VAR                 "4th argument must be function-variable !"
#define RLAB_ERROR_ARG4_MDR                      "4th argument must be real !"
#define RLAB_ERROR_ARG4_MDR_SCALAR               "4th argument must be real scalar !"
#define RLAB_ERROR_ARG4_MDR_VECTOR               "4th argument must be real vector !"
#define RLAB_ERROR_ARG4_MDR_MATRIX               "4th argument must be real matrix !"
#define RLAB_ERROR_ARG4_MDS                      "4th argument must be MDS !"
#define RLAB_ERROR_ARG4_MDS_SCALAR               "4th argument must be string !"
#define RLAB_ERROR_ARG4_MDS_VECTOR               "4th argument must be string vector !"
#define RLAB_ERROR_ARG4_MDS_MATRIX               "4th argument must be string matrix !"

#define RLAB_ERROR_ARG5_FUNC_VAR                 "5th argument must be function-variable !"
#define RLAB_ERROR_ARG5_MDR                      "5th argument must be real !"
#define RLAB_ERROR_ARG5_MDR_SCALAR               "5th argument must be real scalar !"
#define RLAB_ERROR_ARG5_MDR_VECTOR               "5th argument must be real vector !"
#define RLAB_ERROR_ARG5_MDR_MATRIX               "5th argument must be real matrix !"
#define RLAB_ERROR_ARG5_MDS                      "5th argument must be MDS !"
#define RLAB_ERROR_ARG5_MDS_SCALAR               "5th argument must be string !"
#define RLAB_ERROR_ARG5_MDS_VECTOR               "5th argument must be string vector !"
#define RLAB_ERROR_ARG5_MDS_MATRIX               "5th argument must be string matrix !"

#define RLAB_ERROR_ARG6_FUNC_VAR                 "6th argument must be function-variable !"
#define RLAB_ERROR_ARG6_MDR                      "6th argument must be real !"
#define RLAB_ERROR_ARG6_MDR_SCALAR               "6th argument must be real scalar !"
#define RLAB_ERROR_ARG6_MDR_VECTOR               "6th argument must be real vector !"
#define RLAB_ERROR_ARG6_MDR_MATRIX               "6th argument must be real matrix !"
#define RLAB_ERROR_ARG6_MDS                      "6th argument must be MDS !"
#define RLAB_ERROR_ARG6_MDS_SCALAR               "6th argument must be string !"
#define RLAB_ERROR_ARG6_MDS_VECTOR               "6th argument must be string vector !"
#define RLAB_ERROR_ARG6_MDS_MATRIX               "6th argument must be string matrix !"

#define RLAB_ERROR_ARG7_FUNC_VAR                 "7th argument must be a function-variable !"

#define RLAB_ERROR_ARG8_FUNC_VAR                 "8th argument must be a function-variable !"
#define RLAB_ERROR_ARG8_MDS_VECTOR               "8th argument must be string vector !"

#define RLAB_ERROR_ARG9_FUNC_VAR                 "9th argument must be a function-variable !"
#define RLAB_ERROR_ARG9_MDR_MATRIX               "9th argument must be real matrix !"

#define RLAB_ERROR_ARG10_FUNC_VAR                "10th argument must be a function-variable !"

#define RLAB_ERROR_TERRIBLE_INTERNAL_ERROR       "Terrible internal error occured. Cannot continue !"
#define RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR "Jacobian function must return MDR !"
#define RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM   "Incorrectly dimensioned jacobian function !"
#define RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR      "Rhs function must return MDR !"
#define RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR_OR_MDC "Rhs function must return MDR or MDC!"
#define RLAB_ERROR_RHS_FUNC_INCORRECT_DIM        "Incorrectly dimensioned rhs vector !"
#define RLAB_ERROR_BV_FUNC_MUST_RETURN_MDR       "Boundary-value function must return MDR !"
#define RLAB_ERROR_BV_FUNC_INCORRECT_DIM         "Incorrectly dimensioned boundary-value function !"
#define RLAB_ERROR_BV_JAC_INCORRECT_DIM          "Incorrectly dimensioned boundary-value jacobian !"
#define RLAB_ERROR_ONE_ARG_REQUIRED              "One argument required !"
#define RLAB_ERROR_AT_LEAST_ONE_ARG_REQUIRED     "At least one argument required !"
#define RLAB_ERROR_TWO_ARG_REQUIRED              "Two arguments required !"
#define RLAB_ERROR_AT_LEAST_TWO_ARG_REQUIRED     "At least two arguments required !"
#define RLAB_ERROR_AT_LEAST_THREE_ARG_REQUIRED   "At least three arguments required !"
#define RLAB_ERROR_NONE_OR_ONE_ARG_REQUIRED      "None or one argument required !"
#define RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED       "One or two arguments required !"
#define RLAB_ERROR_ONE_OR_THREE_ARG_REQUIRED     "One or three arguments required !"
#define RLAB_ERROR_ONE_OR_FOUR_ARG_REQUIRED      "One or four arguments required !"
#define RLAB_ERROR_ONE_TO_THREE_ARG_REQUIRED     "One, two or three arguments required !"
#define RLAB_ERROR_ONE_TO_FOUR_ARG_REQUIRED      "One to four arguments required !"
#define RLAB_ERROR_ONE_TO_FIVE_ARG_REQUIRED      "One to five arguments required !"
#define RLAB_ERROR_TWO_TO_THREE_ARG_REQUIRED     "Two or three arguments required !"
#define RLAB_ERROR_TWO_TO_FOUR_ARG_REQUIRED      "Two to four arguments required !"
#define RLAB_ERROR_TWO_TO_FIVE_ARG_REQUIRED      "Two to five arguments required !"
#define RLAB_ERROR_TWO_OR_FOUR_ARG_REQUIRED      "Two or four arguments required !"
#define RLAB_ERROR_TWO_3_5_OR_6_ARG_REQUIRED     "Two, three, five or six arguments required !"
#define RLAB_ERROR_THREE_ARG_REQUIRED            "Three arguments required !"
#define RLAB_ERROR_AT_MOST_THREE_ARG_REQUIRED    "At most three arguments required !"
#define RLAB_ERROR_ONE_TO_FOUR_ARG_REQUIRED      "One to four arguments required !"
#define RLAB_ERROR_AT_LEAST_FOUR_ARG_REQUIRED    "At least four arguments required !"
#define RLAB_ERROR_FOUR_ARG_REQUIRED             "Four arguments required !"
#define RLAB_ERROR_ONE_TO_FIVE_ARG_REQUIRED      "One to five arguments required !"
#define RLAB_ERROR_AT_LEAST_FIVE_ARGS_REQUIRED   "At least five arguments required !"
#define RLAB_ERROR_NINE_ARG_REQUIRED             "Nine arguments required !"

#define RLAB_ERROR_ARG_ENTRY_MDR_VECTOR          "Entry must be real vector !"
#define RLAB_ERROR_LSFIT_TOO_FEW_ENTRIES         "Sample size too small for least-square fit !"

#define RLAB_ERROR_SOLVER_NO_RECURSIVE_CALLS     "Solver does not allow recursive calls !"

#define RLAB_ERROR_CANNOT_OPEN_FILE_FOR_READ     "Cannot open file for read !"

#define RLAB_ERROR_ARG1_ARG2_SAME_NUMBER_ROWS    "ARG1 and ARG2 must have same number of rows !"
#define RLAB_ERROR_ARG1_ARG2_SAME_NUMBER_COLS    "ARG1 and ARG2 must have same number of columns !"
#define RLAB_ERROR_ARG1_ARG2_SAME_LENGTH         "ARG1 and ARG2 must have same length !"
#define RLAB_ERROR_ARG1_ARG2_SAME_TYPE           "ARG1 and ARG2 must have same type !"
#define RLAB_ERROR_ARG1_STAT                     "ARG1 must be real vector or list <<val;wgt>> !"
#define RLAB_ERROR_ARG1_MISMATCH_VAL_WGT         "ARG1 list <<val;wgt>> has mismatched entries !"
#define RLAB_ERROR_ARG2_STAT                     "ARG2 must be real vector or list <<val;wgt>> !"
#define RLAB_ERROR_ARG2_MISMATCH_VAL_WGT         "ARG2 list <<val;wgt>> has mismatched entries !"
#define RLAB_ERROR_ARG3_STAT                     "ARG3 must be real vector or list <<val;wgt>> !"
#define RLAB_ERROR_ARG3_MISMATCH_VAL_WGT         "ARG3 list <<val;wgt>> has mismatched entries !"

#define RLAB_ERROR_OUT_OF_MEMORY                  "Out of memory"
#define RLAB_ERROR_FPE                            "Floating Point Exception"
#define RLAB_ERROR_OPERATION_FAILED               "operation not supported"
#define RLAB_ERROR_PRINT_FAILED                   "Print operation not supported"
#define RLAB_ERROR_APPEND_FAILED                  "Append operation not supported"
#define RLAB_ERROR_APPEND_REQUIRES                "Append operation supported only on row-vectors or row-like matrices"
#define RLAB_ERROR_STACK_FAILED                   "Stack operation not supported"
#define RLAB_ERROR_STACK_REQUIRES                 "Stack operation supported only on column-vectors or column-like matrices"
#define RLAB_ERROR_INDEX_OUT_OF_BOUNDS            "Index out of bounds"
#define RLAB_ERROR_SUBMAT_ASSIGNMENT_FAILED       "Sub-matrix assignment operation not supported"
#define RLAB_ERROR_SUBMAT_EVAL_FAILED             "Sub-matrix index evaluation failed"
#define RLAB_ERROR_SUBMATVEC_ASSIGNMENT_FAILED    "Sub-matrix-vector assignment operation not supported"
#define RLAB_ERROR_SUBMAT_INDEX_COERCION_FAILED   "Sub-matrix index coercion operation not supported"
#define RLAB_ERROR_FORLOOPVAL_FAILED              "For-loop-value operation not supported for entity type"
#define RLAB_ERROR_FORLOOP_CANNOT_ITERATE         "Cannot iterate through For-loop"
#define RLAB_ERROR_OPMAT_INVALID_OBJECT_ASSIGN    "Invalid object for vector/matrix-assign"
#define RLAB_ERROR_SUBVEC_INDEX_COERCION_FAILED   "Sub-vector index coercion operation not supported"
#define RLAB_ERROR_SUBVEC_OPERATION_FAILED        "Sub-vector operation not supported"

#define RLAB_WARNING_READM_GENERIC                "Warning: Perhaps file contains no data!"
#define RLAB_WARNING_IGNORE                       "Warning: Incorrect input parameters. Command is being ignored!"
#define RLAB_ERROR_FUNC_UNDEFINED                 "No such function!"
#define RLAB_ERROR_FUNC_CANT_USE_UNDEF_VAR        "Cannot use undefined variable as function!"
#define RLAB_ERROR_FUNC_CANT_USE_CONST            "Cannot use numeric constant as function!"
#define RLAB_ERROR_FUNC_INVALID                   "Invalid function variable!"
#define RLAB_ERROR_HIE                            "Horrible Internal Error"
#define RLAB_ERROR_NO_FUNC_SYM_TABLE              "No Function symbol table!"
//
//
//
//
//
//
//
// embedded languages, et c.
#define RLAB_NAME_EMBED_SPICE   "spice"
#define RLAB_NAME_EMBED_PYTHON  "py"
#define RLAB_NAME_EMBED_JVM     "jvm"
#define RLAB_NAME_EMBED_GPHOTO2 "gp"
#define RLAB_NAME_EMBED_IM      "image"
#define RLAB_NAME_EMBED_GLPK    "lpx"

// rfft:
#define RLAB_NAME_FILTER_Y0       "y0"
#define RLAB_NAME_FILTER_ZI       "zi"
#define RLAB_NAME_FILTER_PERIODIC "periodic"


// generic names for solver parameters:
#define RLAB_NAME_GEN_IMETHOD "imethod"
#define RLAB_NAME_GEN_METHOD  "method"
#define RLAB_NAME_GEN_STDOUT  "stdout"
#define RLAB_NAME_GEN_MAXITER "maxi"
#define RLAB_NAME_GEN_EREL    "erel"
#define RLAB_NAME_GEN_EABS    "eabs"
#define RLAB_NAME_GEN_TOL     "tol"
#define RLAB_NAME_GEN_EPS     "eps"
#define RLAB_NAME_GEN_COEF    "coef"
#define RLAB_NAME_GEN_NAME    "name"
#define RLAB_NAME_GEN_VECTOR  "vec"
#define RLAB_NAME_GEN_INDEX   "idx"
#define RLAB_NAME_GEN_VALUE   "val"
#define RLAB_NAME_GEN_WEIGHT  "wgt"
#define RLAB_NAME_GEN_SIZE    "size"
#define RLAB_NAME_GEN_DEGREE  "deg"
#define RLAB_NAME_GEN_ROW     "row"
#define RLAB_NAME_GEN_COL     "col"
#define RLAB_NAME_GEN_X       "x"
#define RLAB_NAME_GEN_Y       "y"
#define RLAB_NAME_GEN_WIDTH   "width"
#define RLAB_NAME_GEN_BOUNDS  "bounds"
#define RLAB_NAME_GEN_CONSTRAINTS "constraints"



// odeiv
#define RLAB_NAME_ODEIV_STEP    "step"
#define RLAB_NAME_ODEIV_AY      "ay"
#define RLAB_NAME_ODEIV_ADYDT   "adydt"
#define RLAB_NAME_ODEIV_EREL     RLAB_NAME_GEN_EREL
#define RLAB_NAME_ODEIV_EABS     RLAB_NAME_GEN_EABS
#define RLAB_NAME_ODEIV_DELTA_T "delta_t"
#define RLAB_NAME_ODEIV_PHASE   "phase_space"

// odebv
#define RLAB_NAME_ODEBV_NMAX    "x_maxres"
#define RLAB_NAME_ODEBV_TOL      RLAB_NAME_GEN_TOL

// odae
#define RLAB_NAME_ODAE_IMETHOD   RLAB_NAME_GEN_IMETHOD
#define RLAB_NAME_ODAE_STEP     "step"
#define RLAB_NAME_ODAE_INDEX1   "index1"
#define RLAB_NAME_ODAE_INDEX2   "index2"
#define RLAB_NAME_ODAE_INDEX3   "index3"
#define RLAB_NAME_ODAE_EABS      RLAB_NAME_GEN_EABS
#define RLAB_NAME_ODAE_EREL      RLAB_NAME_GEN_EREL
#define RLAB_NAME_ODAE_MAXITER   RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_ODAE_ORDER    "order"
#define RLAB_NAME_ODAE_STDOUT    RLAB_NAME_GEN_STDOUT

// odesl
#define RLAB_NAME_ODESL_TOL      RLAB_NAME_GEN_EABS
#define RLAB_NAME_ODESL_STDOUT   RLAB_NAME_GEN_STDOUT

// pdecol, bacol
#define RLAB_NAME_PDECOL_KORD    "kord"
#define RLAB_NAME_PDECOL_EREL     RLAB_NAME_GEN_EREL
#define RLAB_NAME_PDECOL_DT      "dt"
#define RLAB_NAME_PDECOL_NCC     "ncc"
#define RLAB_NAME_PDECOL_IODE    "iode"
#define RLAB_NAME_PDECOL_MAXDER  "maxder"
#define RLAB_NAME_PDECOL_NOGAUSS "nogauss"
#define RLAB_NAME_PDECOL_STDOUT   RLAB_NAME_GEN_STDOUT
#define RLAB_NAME_BACOL_KORD     "kord"
#define RLAB_NAME_BACOL_IDIR     "idir"
#define RLAB_NAME_BACOL_IMAX     "imax"
#define RLAB_NAME_BACOL_ISOVER   "isover"
#define RLAB_NAME_BACOL_DT       "dt"
#define RLAB_NAME_BACOL_EREL      RLAB_NAME_GEN_EREL
#define RLAB_NAME_BACOL_EABS      RLAB_NAME_GEN_EABS
#define RLAB_NAME_BACOL_STDOUT    RLAB_NAME_GEN_STDOUT

// asa
#define RLAB_NAME_ASA_PT        "param_type"
#define RLAB_NAME_ASA_INITPTEMP "init_param_temp"
#define RLAB_NAME_ASA_CURV      "curvature"
#define RLAB_NAME_ASA_DX        "dx"
#define RLAB_NAME_ASA_SEQ       "seq"
#define RLAB_NAME_ASA_AFM       "acc_freq_mod"
#define RLAB_NAME_ASA_GFM       "gen_freq_mod"
#define RLAB_NAME_ASA_COSTP     "cost_prec"
#define RLAB_NAME_ASA_MAXCOSTR  "max_cost_repeat"
#define RLAB_NAME_ASA_NOCOSTSM  "no_cost_samples"
#define RLAB_NAME_ASA_TEMPRATS  "temp_ratio_scale"
#define RLAB_NAME_ASA_COSTPSRAT "cost_param_scale_ratio"
#define RLAB_NAME_ASA_TEMPANNS  "temp_ann_scale"
#define RLAB_NAME_ASA_LIMITACC  "limit_acc"
#define RLAB_NAME_ASA_LIMITGEN  "limit_gen"
#define RLAB_NAME_ASA_LIMIGS    "limit_inv_gen_states"
#define RLAB_NAME_ASA_AGR       "acc_to_gen_ratio"

// diffevol
#define RLAB_NAME_DIFFEVOL_STRATEGY "strategy"
#define RLAB_NAME_DIFFEVOL_F        "f"
#define RLAB_NAME_DIFFEVOL_FDITHER  "f_dither"
#define RLAB_NAME_DIFFEVOL_CR       "cr"
#define RLAB_NAME_DIFFEVOL_GENMAX   "genmax"
#define RLAB_NAME_DIFFEVOL_IRNG     "irng"
#define RLAB_NAME_DIFFEVOL_STDOUT    RLAB_NAME_GEN_STDOUT
#define RLAB_NAME_DIFFEVOL_IPRINT   "iprint"
#define RLAB_NAME_DIFFEVOL_TARGET   "target"
#define RLAB_NAME_DIFFEVOL_TERR     "terr"
#define RLAB_NAME_DIFFEVOL_MFC      "max_fail_count"
#define RLAB_NAME_DIFFEVOL_EABS      RLAB_NAME_GEN_EABS
#define RLAB_NAME_DIFFEVOL_BOUNDS    RLAB_NAME_GEN_BOUNDS
#define RLAB_NAME_DIFFEVOL_XWIDTH    RLAB_NAME_GEN_WIDTH
#define RLAB_NAME_DIFFEVOL_WIDTHWGT "width_wgt"

// rng / gsl
#define RLAB_NAME_RNG_NAME      "name"
#define RLAB_NAME_RNG_STATE     "state"

// statistics
#define RLAB_NAME_STAT_VALUE     RLAB_NAME_GEN_VALUE
#define RLAB_NAME_STAT_WEIGHT    RLAB_NAME_GEN_WEIGHT
#define RLAB_NAME_STAT_USEIDX     "use_datum"
#define RLAB_NAME_STAT_IGNOREIDX  "ignore_datum"
#define RLAB_NAME_STAT_IGNOREINF  "ignore_inf"
#define RLAB_NAME_STAT_ROWDOM     "row_dominant"
#define RLAB_NAME_STAT_BIAS       "bias"
#define RLAB_NAME_STAT_MEAN       "mean"
#define RLAB_NAME_STAT_FACTOR     "factor"
#define RLAB_NAME_STAT_FLAT       "flat"
#define RLAB_NAME_STAT_METRIC   "metric"
#define RLAB_NAME_STAT_FEATURE  "feature"
#define RLAB_NAME_STAT_DISTMAT  "distmat"
#define RLAB_NAME_HIST1D_BIN    "bin"
#define RLAB_NAME_HIST1D_RANGE  "range"
#define RLAB_NAME_HIST1D_INFO   "_info"
#define RLAB_NAME_HIST1D_PINFS       "pos_inf"
#define RLAB_NAME_HIST1D_NINFS       "neg_inf"
#define RLAB_NAME_HIST1D_NANS        "nan"
#define RLAB_NAME_HIST1D_TRASH_MAX   "sup_max"
#define RLAB_NAME_HIST1D_TRASH_MIN   "sub_min"
#define RLAB_NAME_HIST2D_BIN     RLAB_NAME_HIST1D_BIN
#define RLAB_NAME_HIST2D_XRANGE "xrange"
#define RLAB_NAME_HIST2D_YRANGE "yrange"


// randomize : lapack / gsl
#define RLAB_NAME_RAND_BWUP     "bandwidth_upper"
#define RLAB_NAME_RAND_BWLO     "bandwidth_lower"
#define RLAB_NAME_RAND_SPAR     "sparsity"
#define RLAB_NAME_RAND_SYMM     "symmetric"
#define RLAB_NAME_RAND_HERM     "hermitean"
#define RLAB_NAME_RAND_SVAL     "singular_values"

// siman / gsl
#define RLAB_NAME_GSL_SIMAN_TI  "Ti"
#define RLAB_NAME_GSL_SIMAN_TF  "Tf"
#define RLAB_NAME_GSL_SIMAN_MU  "mu"
#define RLAB_NAME_GSL_SIMAN_NT  "nT"

// monte carlo integration / gsl
#define RLAB_NAME_GSL_MC_NCALL  "ncalls"
#define RLAB_NAME_GSL_MC_NTH    "ntherm"
#define RLAB_NAME_GSL_MC_CHI    "chicomp"
#define RLAB_NAME_GSL_MC_MODE   "mode"
#define RLAB_NAME_GSL_MC_ALPHA  "alpha"
#define RLAB_NAME_GSL_MC_EF     "efrac"
#define RLAB_NAME_GSL_MC_DITH   "dither"
#define RLAB_NAME_GSL_MC_EABS    RLAB_NAME_GEN_EABS

//
// glpk:
//
//  input
#define RLAB_NAME_GLPK_STRUCT_OPTDIR      "opt_direction"
#define RLAB_NAME_GLPK_STRUCT_OBJECTIVE   "objective"
#define RLAB_NAME_GLPK_STRUCT_CONSTRAINT   RLAB_NAME_GEN_CONSTRAINTS
#define RLAB_NAME_GLPK_STRUCT_BOUNDS_AUX  "bounds_row"
#define RLAB_NAME_GLPK_STRUCT_BOUNDS_STR  "bounds_col"
#define RLAB_NAME_GLPK_STRUCT_CLASS       "problem"
#define RLAB_NAME_GLPK_STRUCT_COLS_BIN    "col_bin"
#define RLAB_NAME_GLPK_STRUCT_COLS_INT    "col_int"
#define RLAB_NAME_GLPK_STRUCT_C0          "c0"
// input
#define RLAB_NAME_GLPK_PRINT_SOL  "print_sol"
#define RLAB_NAME_GLPK_METHOD      RLAB_NAME_GEN_METHOD
#define RLAB_NAME_GLPK_ORDER      "ordering"
#define RLAB_NAME_GLPK_OBJMAX     "obj_max"
#define RLAB_NAME_GLPK_OBJMIN     "obj_min"
#define RLAB_NAME_GLPK_SMPRE      "presolve"
#define RLAB_NAME_GLPK_MIPBT      "mip_btrack"
#define RLAB_NAME_GLPK_MAXI        RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_GLPK_PRICE      "pricing"
#define RLAB_NAME_GLPK_RATIO      "r_test"
#define RLAB_NAME_GLPK_TOLBND     "tol_bnd"
#define RLAB_NAME_GLPK_TOLDJ      "tol_dj"
#define RLAB_NAME_GLPK_TOLPIV     "tol_piv"
#define RLAB_NAME_GLPK_STDOUT      RLAB_NAME_GEN_STDOUT
#define RLAB_NAME_GLPK_MIP_TOLINT   "tol_int"
#define RLAB_NAME_GLPK_MIP_TOLOBJ   "tol_obj"
#define RLAB_NAME_GLPK_MIP_BRA      "branch"
#define RLAB_NAME_GLPK_MIP_BTR      "backtrack"
#define RLAB_NAME_GLPK_MIP_PP       "preprocess"
#define RLAB_NAME_GLPK_MIP_FP       "feas_pump"
#define RLAB_NAME_GLPK_MIP_PS       "proximity_search"
#define RLAB_NAME_GLPK_MIP_GMI      "gmi_cuts"
#define RLAB_NAME_GLPK_MIP_MIR      "mir_cuts"
#define RLAB_NAME_GLPK_MIP_COV      "cov_cuts"
#define RLAB_NAME_GLPK_MIP_CLQ      "clq_cuts"
#define RLAB_NAME_GLPK_MIP_MIPGAP   "mip_gap"
#define RLAB_NAME_GLPK_MIP_CB       "cb_size"
#define RLAB_NAME_GLPK_MIP_PRESOLVE "presolve"
#define RLAB_NAME_GLPK_MIP_BINARIZE "binarize"

// findroot / gsl
#define RLAB_NAME_GSL_ROOT1_EABS    RLAB_NAME_GEN_EABS
#define RLAB_NAME_GSL_ROOT1_EREL    RLAB_NAME_GEN_EREL
#define RLAB_NAME_GSL_ROOT1_METHOD  RLAB_NAME_GEN_IMETHOD
#define RLAB_NAME_GSL_ROOT1_MAXI    RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_GSL_ROOT1_RHS     "rhs"

// findroots / gsl, hompack, pbundle
#define RLAB_NAME_GSL_ROOTS_EABS    RLAB_NAME_GEN_EABS
#define RLAB_NAME_GSL_ROOTS_EREL    RLAB_NAME_GEN_EREL
#define RLAB_NAME_GSL_ROOTS_METHOD  RLAB_NAME_GEN_IMETHOD
#define RLAB_NAME_GSL_ROOTS_MAXI    RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_GSL_ROOTS_STDOUT  RLAB_NAME_GEN_STDOUT
// hompack / findroots
#define RLAB_NAME_GSL_ROOTS_ARCRE  "arcre"
#define RLAB_NAME_GSL_ROOTS_ARCAE  "arcae"
#define RLAB_NAME_GSL_ROOTS_ANSRE  "ansre"
#define RLAB_NAME_GSL_ROOTS_ANSAE  "ansae"
#define RLAB_NAME_GSL_ROOTS_FACA   "faca"
#define RLAB_NAME_GSL_ROOTS_CURS   "curs"
#define RLAB_NAME_GSL_ROOTS_IREP   "irep"
#define RLAB_NAME_GSL_ROOTS_STEP   "step"

// leastsquares / gsl
#define RLAB_NAME_GSL_LS_MAXDEG    "maxdeg"
#define RLAB_NAME_GSL_LS_F         "f"
#define RLAB_NAME_GSL_LS_EPS       "eps"
#define RLAB_NAME_GSL_LS_EABS       RLAB_NAME_GEN_EABS
#define RLAB_NAME_GSL_LS_EREL       RLAB_NAME_GEN_EREL
#define RLAB_NAME_GSL_LS_CONVT     "convtest"
#define RLAB_NAME_GSL_LS_MAXITER    RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_GSL_LS_STDOUT     RLAB_NAME_GEN_STDOUT

// odrpack
#define RLAB_NAME_ODR_FIXP         "fix_p"
#define RLAB_NAME_ODR_TAUFAC       "taufac"
#define RLAB_NAME_ODR_SSTOL        "sstol"
#define RLAB_NAME_ODR_PARTOL       "partol"
#define RLAB_NAME_ODR_SCALEP       "scale_p"
#define RLAB_NAME_ODR_IMETHOD       RLAB_NAME_GEN_IMETHOD
#define RLAB_NAME_ODR_MAXITER       RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_ODR_STDOUT        RLAB_NAME_GEN_STDOUT
#define RLAB_NAME_ODR_BOUNDS        RLAB_NAME_GEN_BOUNDS

// min / gsl
#define RLAB_NAME_GSL_MIN1_EABS     RLAB_NAME_GEN_EABS
#define RLAB_NAME_GSL_MIN1_EREL     RLAB_NAME_GEN_EREL
#define RLAB_NAME_GSL_MIN1_MAXITER  RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_GSL_MIN1_IMETHOD  RLAB_NAME_GEN_IMETHOD

// mins / gsl
// methods requiring no derivatives
#define RLAB_NAME_GSL_MINS_EABS     RLAB_NAME_GEN_EABS
#define RLAB_NAME_GSL_MINS_SS      "ss"
#define RLAB_NAME_GSL_MINS_MAXITER  RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_GSL_MINS_STDOUT   RLAB_NAME_GEN_STDOUT
#define RLAB_NAME_GSL_MINS_IMETHOD  RLAB_NAME_GEN_IMETHOD
#define RLAB_NAME_NUO_RHOBEG       "rhobeg"
#define RLAB_NAME_NUO_RHOEND       "rhoend"
#define RLAB_NAME_POW_BOUNDS        RLAB_NAME_GEN_BOUNDS
#define RLAB_NAME_POW_CONSTR        RLAB_NAME_GEN_CONSTRAINTS
#define RLAB_NAME_NUO_NTP          "ntp"
// methods with derivatives
#define RLAB_NAME_PB_MINS_TOLGDSP   RLAB_NAME_GEN_TOL
#define RLAB_NAME_GSL_MINS_STEP    "step"
#define RLAB_NAME_GSL_MINS_TARGET  "target"
#define RLAB_NAME_GSL_MINS_STANDSTILL "standstill"
#define RLAB_NAME_CNM_NRKSTEP      "nrkstep"
#define RLAB_NAME_CNM_TCON         "tolcon"
#define RLAB_NAME_CNM_ISTEP        "istep"
#define RLAB_NAME_CNM_ETA          "eta"
#define RLAB_NAME_CNM_MET          "met"
#define RLAB_NAME_CNM_MES          "mes"
#define RLAB_NAME_CNM_CONVX        "convx"
#define RLAB_NAME_CNM_CONVF        "convf"
#define RLAB_NAME_CNM_IMETHOD       RLAB_NAME_GEN_IMETHOD
#define RLAB_NAME_CNM_MAXITER       RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_CNM_STDOUT        RLAB_NAME_GEN_STDOUT

// nintegrate / gsl, genzpak
#define RLAB_NAME_GSL_INVLT_EREL    RLAB_NAME_GEN_EREL
#define RLAB_NAME_GSL_INVLT_SIGMA  "sigma"
#define RLAB_NAME_GSL_INVLT_SSBAR  "ssbar"
#define RLAB_NAME_GSL_INVLT_FEVAL  "maxfeval"
#define RLAB_NAME_GSL_NINT_EABS     RLAB_NAME_GEN_EABS
#define RLAB_NAME_GSL_NINT_EREL     RLAB_NAME_GEN_EREL
#define RLAB_NAME_GSL_NINT_IKEY    "ikey"
#define RLAB_NAME_GSL_NINT_MAXITER  RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_GSL_NINT_MAXINT  "maxn"
#define RLAB_NAME_GENZ_IMETHOD      RLAB_NAME_GEN_IMETHOD
#define RLAB_NAME_GENZ_IKEY        "ikey"
#define RLAB_NAME_GENZ_MAXI         RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_GENZ_MINI        "mini"
#define RLAB_NAME_GENZ_EABS         RLAB_NAME_GEN_EABS
#define RLAB_NAME_GENZ_EREL         RLAB_NAME_GEN_EREL
#define RLAB_NAME_GSL_CONV_EABS     RLAB_NAME_GEN_EABS
#define RLAB_NAME_GSL_CONV_EREL     RLAB_NAME_GEN_EREL
#define RLAB_NAME_GSL_CONV_IKEY     RLAB_NAME_GSL_NINT_IKEY
#define RLAB_NAME_GSL_CONV_MAXI     RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_GSL_CONV_RES     "res"

// wavelet / gsl
#define RLAB_NAME_GSL_WV_FAMILY    "family"
#define RLAB_NAME_GSL_WV_INDEX     "index"
#define RLAB_NAME_GSL_WV_PROC      "proc_standard"

// harminv
#define RLAB_NAME_HARM_AMPABS      "amp_abs"
#define RLAB_NAME_HARM_AMPREL      "amp_rel"
#define RLAB_NAME_HARM_EABS         RLAB_NAME_GEN_EABS
#define RLAB_NAME_HARM_EREL         RLAB_NAME_GEN_EREL
#define RLAB_NAME_HARM_Q           "Q"
#define RLAB_NAME_HARM_DT          "dt"
#define RLAB_NAME_HARM_DENS        "density"
#define RLAB_NAME_HARM_NF          "nf"
#define RLAB_NAME_HARM_SOLVEM      "solve_mode"

// spline
// #define RLAB_NAME_INTERPL_YPL      "yprime_left"
// #define RLAB_NAME_INTERPL_YPR      "yprime_right"
// #define RLAB_NAME_INTERPL_IDX_STATE "store_state"
#define RLAB_NAME_LINTERP_ORDER    "order"

// gcvspline
#define RLAB_NAME_GCV_DEGREE        RLAB_NAME_GEN_DEGREE
#define RLAB_NAME_GCV_SMOOTH       "smooth"
#define RLAB_NAME_GCV_VAR          "var"
#define RLAB_NAME_GCV_DF           "df"

// dierckx
#define RLAB_NAME_DIERCKX_MAXKNOTS "max_knots"
#define RLAB_NAME_DIERCKX_DEGREE    RLAB_NAME_GEN_DEGREE
#define RLAB_NAME_DIERCKX_RANGE    "range"
#define RLAB_NAME_DIERCKX_TOL       RLAB_NAME_GEN_TOL
#define RLAB_NAME_DIERCKX_MAXITER   RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_DIERCKX_RANKTHR  "rank_threshold"
#define RLAB_NAME_DIERCKX_PER      "periodic"
#define RLAB_NAME_DIERCKX_YP1      "yprime_i"
#define RLAB_NAME_DIERCKX_YP2      "yprime_f"
#define RLAB_NAME_DIERCKX_CONVEX   "convex"

// polyroots, nseries
#define RLAB_NAME_POLY_EPS          RLAB_NAME_GEN_EPS
#define RLAB_NAME_POLY_DMIN        "dmin"
#define RLAB_NAME_POLY_DMAX        "dmax"
#define RLAB_NAME_POLY_MAXITER      RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_NSER_EPS          RLAB_NAME_GEN_EPS
#define RLAB_NAME_NSER_R           "r"

// sprann
#define RLAB_NAME_SPRANN_WH        "WienerHopf"
#define RLAB_NAME_SPRANN_NCYC      "cycles"
#define RLAB_NAME_SPRANN_IMETHOD    RLAB_NAME_GEN_IMETHOD
#define RLAB_NAME_SPRANN_RATE      "rate"
#define RLAB_NAME_SPRANN_MOMENT    "momentum"
#define RLAB_NAME_SPRANN_EABS       RLAB_NAME_GEN_EABS
#define RLAB_NAME_SPRANN_TOL        RLAB_NAME_GEN_TOL
#define RLAB_NAME_SPRANN_STEP      "step"
#define RLAB_NAME_SPRANN_THRES     "threshold"
#define RLAB_NAME_SPRANN_MAXITER    RLAB_NAME_GEN_MAXITER
#define RLAB_NAME_SPRANN_PC        "pc"
#define RLAB_NAME_SPRANN_MSE       "mse"
#define RLAB_NAME_SPRANN_STDOUT     RLAB_NAME_GEN_STDOUT

// rfileio
// serial:
#define RLAB_NAME_SERIAL_DPS       "data_parity_stop"
#define RLAB_NAME_SERIAL_BAUD      "speed"
#define RLAB_NAME_SERIAL_FC        "flow_control"
#define RLAB_NAME_SERIAL_HUPCL     "hupcl"
#define RLAB_NAME_SERIAL_RAW       "raw"
#define RLAB_NAME_SERIAL_DEBUG     "debug"
#define RLAB_NAME_SERIAL_EOL       "eol"
// file:
#define RLAB_NAME_FILE_EOL          RLAB_NAME_SERIAL_EOL
#define RLAB_NAME_FILE_CSP         "csp"
#define RLAB_NAME_FILE_FMT         "format"
#define RLAB_NAME_FILE_BS          "buffer_size"
// socket:
#define RLAB_NAME_SOCKET_MODE      "mode"
#define RLAB_NAME_SOCKET_TIMEOUT   "timeout"

// differentiation and integration
#define RLAB_NAME_DIFFINT_TCH_COEF  RLAB_NAME_GEN_COEF
#define RLAB_NAME_DIFFINT_TCH_INT  "interval"
#define RLAB_NAME_DIFFINT_TCH_DEG  "max_degree"
#define RLAB_NAME_DIFFINT_TCH_NAME  RLAB_NAME_GEN_NAME

// writem parameters
#define RLAB_NAME_WRITEM_SETDTR    "dtr"
#define RLAB_NAME_WRITEM_SETRTS    "rts"
#define RLAB_NAME_WRITEM_FORMAT    "format"
#define RLAB_NAME_WRITEM_COLSEP     RLAB_NAME_FILE_CSP
#define RLAB_NAME_WRITEM_ENDOFLINE  RLAB_NAME_SERIAL_EOL
#define RLAB_NAME_WRITEM_NAN       "nan"
#define RLAB_NAME_WRITEM_INF_POS   "inf_pos"
#define RLAB_NAME_WRITEM_INF_NEG   "inf_neg"

// readm parameters
#define RLAB_NAME_READM_SKIPROWS     "skiprows"
#define RLAB_NAME_READM_CSP           RLAB_NAME_FILE_CSP
#define RLAB_NAME_READM_MIN_LINE_LEN "min_len"
#define RLAB_NAME_READM_COMMENT      "comment"
#define RLAB_NAME_READM_NOTE         "note"
#define RLAB_NAME_READM_USECOLS      "use_cols"
#define RLAB_NAME_READM_USEROWS      "use_rows"
#define RLAB_NAME_READM_JOINROWS     "join_rows"
#define RLAB_NAME_READM_JOINCSP      "join_csp"
#define RLAB_NAME_READM_LSTRIP       "lstrip"
#define RLAB_NAME_READM_START        "start"
#define RLAB_NAME_READM_STOP         "stop"
#define RLAB_NAME_READM_GREP         "grep"

// math
#define RLAB_NAME_MATH_OFFSET        "offset"
#define RLAB_NAME_MATH_BIN           "bin"
#define RLAB_NAME_MATH_BASE          "base"
#define RLAB_NAME_MATH_MESH          "mesh"
#define RLAB_NAME_MATH_CLAMP         "clamp"


// clusters, classifiers and other pattern recognition gugas
#define RLAB_NAME_PATTERN_DATA     "data"
#define RLAB_NAME_PATTERN_FEATURE  "feature"
#define RLAB_NAME_PATTERN_SAMMAN_MAP_GUESS "map"
#define RLAB_NAME_FFNN_NODES       "nodes"
#define RLAB_NAME_FFNN_WEIGHTS      RLAB_NAME_GEN_WEIGHT
#define RLAB_NAME_FFNN_FIX_WEIGHTS  "fix_"RLAB_NAME_GEN_WEIGHT
#define RLAB_NAME_FFNN_BIAS         "bias"
#define RLAB_NAME_FFNN_FIX_BIAS     "fix_"RLAB_NAME_FFNN_BIAS
#define RLAB_NAME_FFNN_ACTIVATION   "act"
#define RLAB_NAME_FFNN_TRANSFER     "transf"
#define RLAB_NAME_FFNN_INIT_RAND    "rand"
#define RLAB_NAME_FFNN_SAMANN_SCALE "scale"
#define RLAB_NAME_FFNN_SAMANN_MEAN  "mean"
#define RLAB_NAME_FFNN_SAMANN_STRESS "stress"
#define RLAB_NAME_FFNN_SAMANN_Y_OUT  "output"
#define RLAB_NAME_FFNN_SAMANN_WGT_RCP_DIST    "rdw"
#define RLAB_NAME_FFNN_NSTEP_ERR_COMP         "nstep_err"
#define RLAB_NAME_FFNN_ERR_DELTA              "delta_err"
#define RLAB_NAME_FFNN_ERR_THRESHOLD          "threshold"
#define RLAB_NAME_FFNN_SAMANN_ISOLVER         RLAB_NAME_GEN_IMETHOD
#define RLAB_NAME_FFNN_SAMANN_RAND_PAIRS      "rand_pairs"
#define RLAB_NAME_PATTERN_BACKPROP_MOMENTUM       "lmomentum"
#define RLAB_NAME_PATTERN_BACKPROP_LEARNING_RATE  "lrate"
#define RLAB_NAME_FFNN_SAMANN_SIMPLEX_S0      RLAB_NAME_GSL_MINS_SS
#define RLAB_NAME_FFNN_SAMANN_SIMPLEX_EABS    RLAB_NAME_GEN_EABS
#define RLAB_NAME_FFNN_SAMANN_NTRIES          "ntries"
#define RLAB_NAME_FFNN_MSE                    "mse"

#define RLAB_ERROR_ARG1_FFNN                  "ARG1 must be list <<"RLAB_NAME_PATTERN_DATA";"RLAB_NAME_PATTERN_FEATURE">> !"

// image magick
#define RLAB_NAME_IMG_QRANGE    "qrange"
#define RLAB_NAME_IMG_ROW        RLAB_NAME_GEN_ROW
#define RLAB_NAME_IMG_COL        RLAB_NAME_GEN_COL
#define RLAB_NAME_IMG_HEIGHT    "height"
#define RLAB_NAME_IMG_WIDTH      RLAB_NAME_GEN_WIDTH
#define RLAB_NAME_IMG_X          RLAB_NAME_GEN_X
#define RLAB_NAME_IMG_Y          RLAB_NAME_GEN_Y
#define RLAB_NAME_IMG_PIXEL_MAP "pixel_map"
#define RLAB_NAME_IMG_PIXEL     "pixel"

// sep
#define RLAB_NAME_SEP_BKG_NOISESTD  "noise_std"
#define RLAB_NAME_SEP_BKG_NOISEVAR  "noise_var"
#define RLAB_NAME_SEP_BKG_NOISEMAT  "noise"
#define RLAB_NAME_SEP_BKG_GAIN      "gain"
#define RLAB_NAME_SEP_BKG_MASK      "mask"
#define RLAB_NAME_SEP_BKG_MASK_THR  "mask_thresh"
#define RLAB_NAME_SEP_BKG_1TILE     "tile"
#define RLAB_NAME_SEP_BKG_FILTER    "filter"
#define RLAB_NAME_SEP_BKG_FILTER_THR  "filter_thresh"
#define RLAB_NAME_SEP_XTR_REL_THR     "rel_thresh"
#define RLAB_NAME_SEP_XTR_ABS_THR     "abs_thresh"
#define RLAB_NAME_SEP_XTR_MINAREA     "min_area"
#define RLAB_NAME_SEP_XTR_MAXDEBAREA    "max_area"
#define RLAB_NAME_SEP_XTR_DEBTHRESHTYPE "deblend_type"
#define RLAB_NAME_SEP_XTR_CONV          "conv"
#define RLAB_NAME_SEP_XTR_DEBLEND       "deblend"
#define RLAB_NAME_SEP_XTR_CONTRAST      "contrast"
#define RLAB_NAME_SEP_XTR_CLEAN       "clean"
#define RLAB_NAME_SEP_XTR_PIXEL_STACK "pixel_stack"
#define RLAB_NAME_SEP_XTR_CENTER      "pos"
#define RLAB_NAME_SEP_XTR_SECMOMNT    "m2"
#define RLAB_NAME_SEP_XTR_ELLIPSE     "ellipse"
#define RLAB_NAME_SEP_XTR_FLUX        "flux"
#define RLAB_NAME_SEP_XTR_CFLUX       "cflux"
#define RLAB_NAME_SEP_XTR_PEAK        "peak"
#define RLAB_NAME_SEP_XTR_PEAK_POS    "peak_pos"
#define RLAB_NAME_SEP_XTR_CPEAK       "cpeak"
#define RLAB_NAME_SEP_XTR_CPEAK_POS   "cpeak_pos"
#define RLAB_NAME_SEP_XTR_CENTER_C    "center_conv"
#define RLAB_NAME_SEP_XTR_PIXELS      "pixel_pos"
#define RLAB_NAME_SEP_XTR_SORT_BY     "sort_by"

// gnuplot
#define RLAB_NAME_GNUPLOT_GNUWINS_WINDOWS "win"
#define RLAB_NAME_GNUPLOT_GNUWINS_DEVICE  "dev"
#define RLAB_NAME_GNUPLOT_GNUWINS_ACTIVE  "act"


#endif /* RLAB_PARAMETER_NAMES_ETC */
















