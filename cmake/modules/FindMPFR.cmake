set(MPFR_FIND_REQUIRED)

find_path(MPFR_INCLUDE_DIR mpfr.h
	PATHS $ENV{GMPDIR} $ENV{MPFRDIR} ${INCLUDE_INSTALL_DIR}
)

find_library(MPFR_LIBRARIES mpfr
	PATHS $ENV{GMPDIR} $ENV{MPFRDIR} ${INCLUDE_INSTALL_DIR}	
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPFR DEFAULT_MSG
                                  MPFR_INCLUDE_DIR MPFR_LIBRARIES)