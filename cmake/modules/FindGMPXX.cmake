find_path(GMPXX_INCLUDE_DIR 
                    NAMES gmpxx.h 
                    PATHS $ENV{GMPXXDIR} $ENV{GMPXX_HOME} $ENV{GMPXX_INCLUDE} $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMP_INCLUDE} ${INCLUDE_INSTALL_DIR}
                    )

find_library(GMPXX_LIBRARIES gmpxx
                    PATHS $ENV{GMPXXDIR} $ENV{GMPXX_HOME} $ENV{GMPXX_LIB} $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMP_LIB} ${LIB_INSTALL_DIR}
                    )


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMPXX DEFAULT_MSG
                                    GMPXX_INCLUDE_DIR GMPXX_LIBRARIES)