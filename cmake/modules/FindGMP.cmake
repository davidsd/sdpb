find_path(GMP_INCLUDE_DIR 
                    NAMES gmp.h 
                    PATHS $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMP_INCLUDE} ${INCLUDE_INSTALL_DIR}
                    )

find_library(GMP_LIBRARIES gmp 
                    PATHS $ENV{GMPDIR} $ENV{GMP_HOME} $ENV{GMP_LIB} ${LIB_INSTALL_DIR}
                    )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP DEFAULT_MSG
                                    GMP_INCLUDE_DIR GMP_LIBRARIES)