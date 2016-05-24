# ==============================================================================
# Composite CMake File
# ==============================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
PROJECT(GP_GENERAL)

# ==============================================================================
# project setup ----------------------------------------------------------------
# ==============================================================================

# what languages are we going to use -------------------------------------------
ENABLE_LANGUAGE(Fortran)

# this is a composite project ----------------------------------------------------
SET(COMPOSITE_PROJECT TRUE)

# shared or static build
IF(NOT LIB_BUILD_TYPE)
	SET(LIB_BUILD_TYPE SHARED)
	MESSAGE("==> LIB_BUILD_TYPE is not defined, default (SHARED) is applied.")
ELSEIF("${LIB_BUILD_TYPE}" MATCHES "SHARED" OR
       "${LIB_BUILD_TYPE}" MATCHES "Shared" OR
       "${LIB_BUILD_TYPE}" MATCHES "shared")
	SET(LIB_BUILD_TYPE SHARED)
	MESSAGE("==> LIB_BUILD_TYPE is set to SHARED.")
ELSEIF("${LIB_BUILD_TYPE}" MATCHES "STATIC" OR
       "${LIB_BUILD_TYPE}" MATCHES "Static" OR
       "${LIB_BUILD_TYPE}" MATCHES "static")
	SET(LIB_BUILD_TYPE STATIC)
	MESSAGE("==> LIB_BUILD_TYPE is set to STATIC.")
ELSE()
	MESSAGE(FATAL_ERROR "ERROR: LIB_BUILD_TYPE must be either SHARED or STATIC!")
ENDIF()

# default build type -----------------------------------------------------------
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE RELEASE)
	MESSAGE("==> CMAKE_BUILD_TYPE is not defined, default (RELEASE) is applied.")
ELSEIF("${CMAKE_BUILD_TYPE}" MATCHES "RELEASE" OR 
       "${CMAKE_BUILD_TYPE}" MATCHES "Release" OR
       "${CMAKE_BUILD_TYPE}" MATCHES "release")
	SET(CMAKE_BUILD_TYPE RELEASE)
	MESSAGE("==> CMAKE_BUILD_TYPE is set to RELEASE.")
ELSEIF("${CMAKE_BUILD_TYPE}" MATCHES "DEBUG" OR
       "${CMAKE_BUILD_TYPE}" MATCHES "Debug" OR
       "${CMAKE_BUILD_TYPE}" MATCHES "debug")
	SET(CMAKE_BUILD_TYPE DEBUG)
	MESSAGE("==> CMAKE_BUILD_TYPE is set to DEBUG.")
ELSE()
	MESSAGE(FATAL_ERROR "ERROR: CMAKE_BUILD_TYPE must be either RELEASE or DEBUG!")
ENDIF()

# default precision ------------------------------------------------------------
IF(NOT PREC)
	SET(PREC DOUBLE)
	MESSAGE("==> PREC is not defined, default (DOUBLE) is applied.")
ELSEIF("${PREC}" MATCHES "DOUBLE" OR
       "${PREC}" MATCHES "Double" OR
       "${PREC}" MATCHES "double" OR
       "${PREC}" MATCHES "DP" OR
       "${PREC}" MATCHES "dp")
	ADD_DEFINITIONS(-DDP)
	SET(PREC DOUBLE)
	MESSAGE("==> PREC is set to DOUBLE.")
ELSEIF("${PREC}" MATCHES "QUAD" OR
       "${PREC}" MATCHES "Quad" OR
       "${PREC}" MATCHES "quad" OR
       "${PREC}" MATCHES "QP" OR
       "${PREC}" MATCHES "qp")
	ADD_DEFINITIONS(-DQP)
	SET(PREC QUAD)
	MESSAGE("==> PREC is set to QUAD.")
ELSE()
	MESSAGE(FATAL_ERROR "ERROR: PREC must be either DOUBLE or QUAD!")
ENDIF()

# find BLAS/LAPACK libraries ---------------------------------------------------
IF("${PREC}" MATCHES "DOUBLE")
	IF(NOT BLAS_VENDOR)
		IF(EXISTS $ENV{MKLROOT})
			SET(BLAS_VENDOR "mkl")
			SET(BLAS_LAPACK_LIBRARIES $ENV{MKLROOT})
			MESSAGE("==> MKL found (MKLROOT): ${BLAS_LAPACK_LIBRARIES}")
		ELSE()
			IF(NOT BLAS_LAPACK_LIBRARIES)
				SET(BLAS_VENDOR "intrinsic")
				MESSAGE("==> no BLAS/LAPACK libraries have been found!")
				MESSAGE("==> intrinsic double precision library is used!")
			ELSE()
				SET(BLAS_VENDOR "lib")
			ENDIF()
		ENDIF()
	ELSE()
		IF("${BLAS_VENDOR}" MATCHES "MKL" OR
                   "${BLAS_VENDOR}" MATCHES "mkl" OR
		   "${BLAS_VENDOR}" MATCHES "Intel" OR
		   "${BLAS_VENDOR}" MATCHES "intel")
			SET(BLAS_VENDOR "mkl")
        		IF(EXISTS $ENV{MKLROOT})
				SET(BLAS_LAPACK_LIBRARIES $ENV{MKLROOT})
                		MESSAGE("-- MKL found (MKLROOT): ${BLAS_LAPACK_LIBRARIES}")
			ELSE()
				MESSAGE(FATAL_ERROR "ERROR: MKLROOT is not defined!")
			ENDIF()
		ELSEIF("${BLAS_VENDOR}" MATCHES "ACC" OR
		       "${BLAS_VENDOR}" MATCHES "acc" OR
		       "${BLAS_VENDOR}" MATCHES "Acc")
			INCLUDE(FindBLAS)
			INCLUDE(FindLAPACK)
			IF(${BLAS_FOUND} AND ${LAPACK_FOUND})
				SET(BLAS_VENDOR "acc")
			ELSE()
				SET(BLAS_VENDOR "intrinsic")
                                MESSAGE("==> no BLAS/LAPACK libraries have been found!")
                                MESSAGE("==> intrinsic double precision library is used!")
			ENDIF()
		ELSEIF("${BLAS_VENDOR}" MATCHES "LIB" OR
		       "${BLAS_VENDOR}" MATCHES "lib" OR
		       "${BLAS_VENDOR}" MATCHES "Lib")
			IF(NOT BLAS_LAPACK_LIBRARIES)
				INCLUDE(FindBLAS)
                        	INCLUDE(FindLAPACK)
				IF(${BLAS_FOUND} AND ${LAPACK_FOUND})
					SET(BLAS_VENDOR "lib")
				ELSE()
					SET(BLAS_VENDOR "intrinsic")
                                        MESSAGE("==> no BLAS/LAPACK libraries have been found!")
                                        MESSAGE("==> intrinsic double precision library is used!")
				ENDIF()
			ENDIF()
		ENDIF()
        ENDIF()
ELSE()
	SET(BLAS_VENDOR "intrinsic")
        MESSAGE("==> intrinsic quad precision library is used!")
ENDIF()

IF("${BLAS_VENDOR}" MATCHES "intrinsic")
ELSE()
	ADD_DEFINITIONS(-DEXTERNAL_BLAS_LAPACK)
ENDIF()

# compiler specifications ------------------------------------------------------
IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
	# Intel
	IF("${BLAS_VENDOR}" MATCHES "mkl")
		SET(BLAS_LAPACK_FLAG "-mkl")
	ELSEIF("${BLAS_VENDOR}" MATCHES "lib")
		IF(NOT BLAS_LAPACK_LIBRARIES)
			SET(BLAS_LAPACK_FLAG "-llapack -lblas")
		ELSE()
			SET(BLAS_LAPACK_FLAG "-L${BLAS_LAPACK_LIBRARIES} -llapack -lblas")
		ENDIF()
	ENDIF()
	IF("${LIB_BUILD_TYPE}" MATCHES "SHARED")
 		SET(CMAKE_Fortran_FLAGS_RELEASE "-fpp -traceback -parallel ${BLAS_LAPACK_FLAG} -openmp -O3")
 	        SET(CMAKE_Fortran_FLAGS_DEBUG "-fpp -O0 -g -traceback -check all -warn all ${BLAS_LAPACK_FLAG} -openmp")
	ELSEIF("${LIB_BUILD_TYPE}" MATCHES "STATIC")
 		SET(CMAKE_Fortran_FLAGS_RELEASE "-fpp -traceback -parallel ${BLAS_LAPACK_FLAG} -openmp -fast")
 		SET(CMAKE_Fortran_FLAGS_DEBUG "-fpp -O0 -g -traceback -check all -warn all ${BLAS_LAPACK_FLAG} -openmp")
	ENDIF()
	# Intel library archiver
	FIND_PROGRAM(XIAR xiar)
	IF(XIAR)
		SET(CMAKE_AR "${XIAR}")
	ENDIF(XIAR)
	MARK_AS_ADVANCED(XIAR)
	# Intel linker
	FIND_PROGRAM(XILD xild)
	IF(XILD)
		SET(CMAKE_LINKER "${XILD}")
	ENDIF(XILD)
	MARK_AS_ADVANCED(XILD)
ELSEIF(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
	# GNU
	IF("${BLAS_VENDOR}" MATCHES "mkl")
		SET(BLAS_LAPACK_FLAG "-Wl,--no-as-needed -L${BLAS_LAPACK_LIBRARIES}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm")
	ELSEIF("${BLAS_VENDOR}" MATCHES "lib")
		IF(NOT BLAS_LAPACK_LIBRARIES)
			SET(BLAS_LAPACK_FLAG "-llapack -lblas")
		ELSE()
			SET(BLAS_LAPACK_FLAG "-L${BLAS_LAPACK_LIBRARIES} -llapack -lblas")
		ENDIF()
	ELSEIF("${BLAS_VENDOR}" MATCHES "acc")
		SET(BLAS_LAPACK_FLAG "-framework Accelerate")
	ENDIF()
 	SET(CMAKE_Fortran_FLAGS_RELEASE "-ffree-form -cpp -O3 -fbacktrace -fopenmp ${BLAS_LAPACK_FLAG}")
 	SET(CMAKE_FOTRTAN_FLAGS_DEBUG   "-ffree-form -cpp -O0 -fbacktrace -Wall -Wextra -fopenmp ${BLAS_LAPACK_FLAG}")
ELSE()
	MESSAGE(FATAL_ERROR "ERROR: ${CMAKE_Fortran_COMPILER_ID} is not supported Fortran compiler!")
ENDIF()

IF("${CMAKE_BUILD_TYPE}" MATCHES "RELEASE")
	MESSAGE("==> Flags: " ${CMAKE_Fortran_FLAGS_RELEASE})
ELSEIF("${CMAKE_BUILD_TYPE}" MATCHES "DEBUG")
	MESSAGE("==>Flags: " ${CMAKE_Fortran_FLAGS_DEBUG})
ELSE()
	MESSAGE(FATAL_ERROR "ERROR: CMAKE_BUILD_TYPE must be either RELEASE or DEBUG!")
ENDIF()

MESSAGE("${BLAS_VENDOR}")
# ==============================================================================
# output directories -----------------------------------------------------------
# ==============================================================================

SET(BIN "${CMAKE_BINARY_DIR}/bin")
SET(LIB "${CMAKE_BINARY_DIR}/lib")
SET(MOD "${CMAKE_BINARY_DIR}/mod")
SET(SRC "${CMAKE_BINARY_DIR}/src")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN})
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB})
SET(CMAKE_Fortran_MODULE_DIRECTORY ${MOD})

# ==============================================================================
# prerequisites ----------------------------------------------------------------
# ==============================================================================

SET(PREC_DIR "${SRC}/prerequisites/prec")
SET(PREC_SRC
	${PREC_DIR}/prec.f90
	)
ADD_LIBRARY(prec_lib ${LIB_BUILD_TYPE} ${PREC_SRC})
SET_TARGET_PROPERTIES(prec_lib PROPERTIES OUTPUT_NAME prec)
SET(PREC_LIB_NAME prec_lib)

IF("${BLAS_VENDOR}" MATCHES "intrinsic")
	SET(BLAS_LAPACK_DIR "${SRC}/prerequisites/blas_lapack")
	SET(BLAS_LAPACK_SRC
		${BLAS_LAPACK_DIR}/blas_lapack.f90
	)
	ADD_LIBRARY(blas_lapack_lib ${LIB_BUILD_TYPE} ${BLAS_LAPACK_SRC})
	TARGET_LINK_LIBRARIES(blas_lapack_lib prec_lib)
	SET_TARGET_PROPERTIES(blas_lapack_lib PROPERTIES OUTPUT_NAME blas_lapack)
	SET(BLAS_LAPACK_LIB_NAME blas_lapack_lib)
ENDIF()

SET(CONSTANTS_DIR "${SRC}/prerequisites/constants")
SET(CONSTANTS_SRC
	${CONSTANTS_DIR}/constants.f90
	)

ADD_LIBRARY(constants_lib ${LIB_BUILD_TYPE} ${CONSTANTS_SRC})
TARGET_LINK_LIBRARIES(constants_lib prec_lib)
SET_TARGET_PROPERTIES(constants_lib PROPERTIES OUTPUT_NAME constants)
SET(CONSTANTS_LIB_NAME constants_lib)


SET(VERBOSITY_DIR "${SRC}/prerequisites/verbosity")
SET(VERBOSITY_SRC
	${VERBOSITY_DIR}/verbosity.f90
	)

ADD_LIBRARY(verbosity_lib ${LIB_BUILD_TYPE} ${VERBOSITY_SRC})
TARGET_LINK_LIBRARIES(verbosity_lib prec_lib)
SET_TARGET_PROPERTIES(verbosity_lib PROPERTIES OUTPUT_NAME verbosity)
SET(VERBOSITY_LIB_NAME verbosity_lib)


SET(POSITIONS_DIR "${SRC}/prerequisites/positions")
SET(POSITIONS_SRC
	${POSITIONS_DIR}/positions.f90
	)

ADD_LIBRARY(positions_lib ${LIB_BUILD_TYPE} ${POSITIONS_SRC})
TARGET_LINK_LIBRARIES(positions_lib prec_lib)
SET_TARGET_PROPERTIES(positions_lib PROPERTIES OUTPUT_NAME positions)
SET(POSITIONS_LIB_NAME positions_lib)

SET(CLUSTERING_DIR "${SRC}/prerequisites/clustering")
SET(CLUSTERING_SRC
	${CLUSTERING_DIR}/clustering.f90
	)

ADD_LIBRARY(clustering_lib ${LIB_BUILD_TYPE} ${CLUSTERING_SRC})
TARGET_LINK_LIBRARIES(clustering_lib positions_lib prec_lib)
SET_TARGET_PROPERTIES(clustering_lib PROPERTIES OUTPUT_NAME clustering)
SET(CLUSTERING_LIB_NAME clustering_lib)

SET(OPTIM_DIR "${SRC}/prerequisites/optim")
SET(OPTIM_SRC
	${OPTIM_DIR}/lbfgs.f90
	${OPTIM_DIR}/lbfgsb.f90
	)

ADD_LIBRARY(optim_lib ${LIB_BUILD_TYPE} ${OPTIM_SRC})
TARGET_LINK_LIBRARIES(optim_lib prec_lib)
IF("${BLAS_VENDOR}" MATCHES "intrinsic")
	TARGET_LINK_LIBRARIES(optim_lib blas_lapack_lib)
ENDIF()
SET_TARGET_PROPERTIES(optim_lib PROPERTIES OUTPUT_NAME optim)
SET(OPTIM_LIB_NAME optim_lib)

SET(FPARSER_DIR "${SRC}/prerequisites/fparser")
SET(FPARSER_SRC
	${FPARSER_DIR}/fparser_dat.f90
	${FPARSER_DIR}/fparser.f90
	)

ADD_LIBRARY(fparser_lib ${LIB_BUILD_TYPE} ${FPARSER_SRC})
TARGET_LINK_LIBRARIES(fparser_lib constants_lib prec_lib)
SET_TARGET_PROPERTIES(fparser_lib PROPERTIES OUTPUT_NAME fparser)
SET(FPARSER_LIB_NAME fparser_lib)

SET(GP_BASIC_DIR "${SRC}/prerequisites/gp_basic")
SET(GP_BASIC_SRC
	${GP_BASIC_DIR}/gp_basic_dat.f90
	${GP_BASIC_DIR}/gp_basic.f90
	)

ADD_LIBRARY(gp_basic_lib ${LIB_BUILD_TYPE} ${GP_BASIC_SRC})
TARGET_LINK_LIBRARIES(gp_basic_lib fparser_lib constants_lib prec_lib)
IF("${BLAS_VENDOR}" MATCHES "intrinsic")
	TARGET_LINK_LIBRARIES(gp_basic_lib blas_lapack_lib)
ENDIF()
SET_TARGET_PROPERTIES(gp_basic_lib PROPERTIES OUTPUT_NAME gp_basic)
SET(GP_BASIC_LIB_NAME gp_basic_lib)

# projects ---------------------------------------------------------------------

SET(GP_GENERAL_DIR "${SRC}/projects/gp_general")
SET(GP_GENERAL_SRC
	${GP_GENERAL_DIR}/gp_general_dat.f90
	${GP_GENERAL_DIR}/gp_general_read.f90
	${GP_GENERAL_DIR}/gp_general_main.f90
	${GP_GENERAL_DIR}/gp_general_write.f90
	${GP_GENERAL_DIR}/gp_general.f90
	)

# ==============================================================================
# create executables -----------------------------------------------------------
# ==============================================================================

# define extensions ------------------------------------------------------------
IF("${CMAKE_BUILD_TYPE}" MATCHES "RELEASE")
	SET(EXT_BUILD "")
ELSEIF("${CMAKE_BUILD_TYPE}" MATCHES "DEBUG")
	SET(EXT_BUILD "_debug")
ENDIF()

IF("${PREC}" MATCHES "DOUBLE")
	SET(EXT_PREC "")
ELSEIF("${PREC}" MATCHES "QUAD")
	SET(EXT_PREC "_quad")
ENDIF()

# add executables --------------------------------------------------------------
ADD_EXECUTABLE("gp_general${EXT_BUILD}${EXT_PREC}" ${GP_GENERAL_SRC})
TARGET_LINK_LIBRARIES("gp_general${EXT_BUILD}${EXT_PREC}"
	${PREC_LIB_NAME}
	${BLAS_LAPACK_LIB_NAME}
	${CONSTANTS_LIB_NAME}
	${VERBOSITY_LIB_NAME}
	${POSITIONS_LIB_NAME}
	${CLUSTERING_LIB_NAME}
	${OPTIM_LIB_NAME}
	${FPARSER_LIB_NAME}
	${GP_BASIC_LIB_NAME}
	)
