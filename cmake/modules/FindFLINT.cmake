set(FLINT_FIND_REQUIRED)

find_path(FLINT_INCLUDE_DIR flint/flint.h flint/fmpz.h flint/nmod.h flint/nmod_vec.h flint/fmpz_mat.h
	PATHS $ENV{CPATH} $ENV{GMPDIR} $ENV{MPFRDIR} ${INCLUDE_INSTALL_DIR}
)

find_library(FLINT_LIBRARIES flint
	PATHS $ENV{LIBRARY_PATH}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT DEFAULT_MSG
                                  FLINT_INCLUDE_DIR FLINT_LIBRARIES)