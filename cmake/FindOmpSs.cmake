# Detect installation of OmpSs
# 1) use OMPSS_HOME env var (or)
# 2) Search for mcc/mcxx in PATH

INCLUDE (FindPackageHandleStandardArgs)
MESSAGE (STATUS "Finding OmpSs location ...")

IF (NOT OMPSS_ROOT AND NOT $ENV{OMPSS_HOME} STREQUAL "")
  MESSAGE (STATUS "Found OMPSS_HOME=$ENV{OMPSS_HOME}")
  SET (OMPSS_ROOT "$ENV{OMPSS_HOME}")
  MESSAGE ("OMPSS_ROOT=${OMPSS_ROOT}")
  FIND_FILE (tapenade "${OMPSS_ROOT}/bin" NO_DEFAULT_PATH)
  SET (OMPSS_BIN "${OMPSS_ROOT}/bin")
  IF (NOT OMPSS_BIN)
    MESSAGE (WARNING "Cannot locate mcc/mcxx executable in OMPSS_HOME=$ENV{OMPSS_HOME}")
    SET (OMPSS_FOUND FALSE)
  ELSE ()
    SET (OMPSS_CXX "${OMPSS_BIN}/mcxx")
    SET (OMPSS_C "${OMPSS_BIN}/mcc")
    SET (OMPSS_FOUND TRUE)
  ENDIF ()
ENDIF ()

# search standard paths
IF (NOT OMPSS_ROOT)
  FIND_PATH (OMPSS_BIN NAMES mcc)
  IF (NOT OMPSS_BIN)
    MESSAGE (WARNING "Cannot locate mcc/mcxx executable in path")
    SET (OMPSS_FOUND FALSE)
  ELSE ()
    MESSAGE (STATUS "Found mcc executable in ${OMPSS_BIN}")
    SET (OMPSS_C "${OMPSS_BIN}/mcc")
    SET (OMPSS_CXX "${OMPSS_BIN}/mcxx")
    SET (OMPSS_FOUND TRUE)
  ENDIF ()
ENDIF ()

MARK_AS_ADVANCED (OMPSS_C OMPSS_CXX)
