#
# This file is part of ACADO Toolkit.
#
# ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
# Copyright (C) 2008-2011 by Boris Houska and Hans Joachim Ferreau.
# All rights reserved.
#
# ACADO Toolkit is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# ACADO Toolkit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with ACADO Toolkit; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

################################################################################
#
# Description:
#	ACADO Toolkit package configuration file
#
#	Defines:
#		- Variable: ACADO_INCLUDE_DIRS
#		- Variable: ACADO_LIBRARY_DIRS
#		- Variable: ACADO_STATIC_LIBRARIES
#		- Variable: ACADO_SHARED_LIBRARIES
#
# Authors:
#	Milan Vukov, milan.vukov@esat.kuleuven.be
#
# Year:
#	2011.
#
# NOTE:
#	- This script is for Linux/Unix use only.
#	- Use this script only (and only :)) if you do not want to install ACADO
#		toolkit. If you install ACADO toolkit you do not need this script.
#
#	- PREREQUISITE: sourced acado_env.sh in your ~/.bashrc file. This script
#		fill try to find ACADO folders, libraries etc., but looking for them
#		in enviromental variables.
#
# Usage:
#	- Linux/Unix: TODO
#
################################################################################

################################################################################
#
# Search for package components
#
################################################################################

MESSAGE( STATUS "********************************************************************************" )
MESSAGE( STATUS "Looking for ACADO toolkit: \n" )

#
# Include folders
#
MESSAGE( STATUS "Looking for ACADO toolkit include directories" )
SET( ACADO_INCLUDE_DIRS $ENV{ACADO_ENV_INCLUDE_DIRS} )
IF( ACADO_INCLUDE_DIRS )
	MESSAGE( STATUS "Found ACADO toolkit include directories: ${ACADO_INCLUDE_DIRS} \n" )
	SET( ACADO_INCLUDE_DIRS_FOUND TRUE )
ELSE( ACADO_INCLUDE_DIRS )
	MESSAGE( STATUS "Could not find ACADO toolkit include directories \n" )
ENDIF( ACADO_INCLUDE_DIRS )

#
# Library folders
#
MESSAGE( STATUS "Looking for ACADO toolkit library directories" )
SET( ACADO_LIBRARY_DIRS $ENV{ACADO_ENV_LIBRARY_DIRS} )
IF( ACADO_LIBRARY_DIRS )
	MESSAGE( STATUS "Found ACADO toolkit library directories: ${ACADO_LIBRARY_DIRS} \n" )
	SET( ACADO_LIBRARY_DIRS_FOUND TRUE )
ELSE( ACADO_LIBRARY_DIRS )
	MESSAGE( STATUS "Could not find ACADO toolkit library directories \n" )
ENDIF( ACADO_LIBRARY_DIRS )

#
# Static libs
#
MESSAGE( STATUS "Looking for ACADO toolkit static libraries" )
SET( ACADO_STATIC_LIBS_FOUND TRUE )

FIND_LIBRARY( ACADO_TOOLKIT_LIB
	NAMES acado_toolkit
	PATHS ${ACADO_LIBRARY_DIRS}
	NO_DEFAULT_PATH
)
IF( ACADO_TOOLKIT_LIB )
	MESSAGE( STATUS "Found ACADO static library: acado_toolkit" )
ELSE( ACADO_TOOLKIT_LIB )
	MESSAGE( STATUS "Could not find ACADO static library: acado_toolkit" )
	SET( ACADO_STATIC_LIBS_FOUND FALSE )
ENDIF( ACADO_TOOLKIT_LIB )
		
FIND_LIBRARY( ACADO_QPOASES_LIB
	NAMES acado_qpOASESextras
	PATHS ${ACADO_LIBRARY_DIRS}
	NO_DEFAULT_PATH
)
IF( ACADO_QPOASES_LIB )
	MESSAGE( STATUS "Found ACADO static library: acado_qpOASESextras" )
ELSE( ACADO_QPOASES_LIB )
	MESSAGE( STATUS "Could not find ACADO static library: acado_qpOASESextras" )
	SET( ACADO_STATIC_LIBS_FOUND FALSE )
ENDIF( ACADO_QPOASES_LIB )
		
FIND_LIBRARY( ACADO_CSPARSE_LIB
	NAMES acado_csparse
	PATHS ${ACADO_LIBRARY_DIRS}
	NO_DEFAULT_PATH
)
IF( ACADO_CSPARSE_LIB )
	MESSAGE( STATUS "Found ACADO static library: acado_csparse" )
ELSE( ACADO_CSPARSE_LIB )
	MESSAGE( STATUS "Could not find ACADO static library: acado_csparse" )
	SET( ACADO_STATIC_LIBS_FOUND FALSE )
ENDIF( ACADO_CSPARSE_LIB )
		
SET( ACADO_STATIC_LIBRARIES
	${ACADO_TOOLKIT_LIB} ${ACADO_QPOASES_LIB} ${ACADO_CSPARSE_LIB}
)

IF( ACADO_STATIC_LIBS_FOUND )
	MESSAGE( STATUS "Found ACADO toolkit static libraries\n" )
ELSE()
	MESSAGE( STATUS "Could not find ACADO toolkit static libraries\n" )
ENDIF()

IF( VERBOSE )
	MESSAGE( STATUS "${ACADO_STATIC_LIBRARIES}\n" )
ENDIF()

#
# Shared libs
#
MESSAGE( STATUS "Looking for ACADO toolkit shared libraries" )
SET( ACADO_SHARED_LIBS_FOUND TRUE )

FIND_LIBRARY( ACADO_TOOLKIT_LIB_S
	NAMES acado_toolkit_s
	PATHS ${ACADO_LIBRARY_DIRS}
	NO_DEFAULT_PATH
)
IF( ACADO_TOOLKIT_LIB_S )
	MESSAGE( STATUS "Found ACADO shared library: acado_toolkit_s" )
ELSE( ACADO_TOOLKIT_LIB_S )
	MESSAGE( STATUS "Could not find ACADO shared library: acado_toolkit_s" )
	SET( ACADO_SHARED_LIBS_FOUND FALSE )
ENDIF( ACADO_TOOLKIT_LIB_S )
		
FIND_LIBRARY( ACADO_QPOASES_LIB_S
	NAMES acado_qpOASESextras_s
	PATHS ${ACADO_LIBRARY_DIRS}
	NO_DEFAULT_PATH
)
IF( ACADO_QPOASES_LIB_S )
	MESSAGE( STATUS "Found ACADO shared library: acado_qpOASESextras_s" )
ELSE( ACADO_QPOASES_LIB_S )
	MESSAGE( STATUS "Could not find ACADO shared library: acado_qpOASESextras_s" )
	SET( ACADO_SHARED_LIBS_FOUND FALSE )
ENDIF( ACADO_QPOASES_LIB_S )
		
FIND_LIBRARY( ACADO_CSPARSE_LIB_S
	NAMES acado_csparse_s
	PATHS ${ACADO_LIBRARY_DIRS}
	NO_DEFAULT_PATH
)
IF( ACADO_CSPARSE_LIB_S )
	MESSAGE( STATUS "Found ACADO shared library: acado_csparse_s" )
ELSE( ACADO_CSPARSE_LIB_S )
	MESSAGE( STATUS "Could not find ACADO shared library: acado_csparse_s" )
	SET( ACADO_SHARED_LIBS_FOUND FALSE )
ENDIF( ACADO_CSPARSE_LIB_S )
		
SET( ACADO_SHARED_LIBRARIES
	${ACADO_TOOLKIT_LIB_S} ${ACADO_QPOASES_LIB_S} ${ACADO_CSPARSE_LIB_S}
)

IF( ACADO_SHARED_LIBS_FOUND )
	MESSAGE( STATUS "Found ACADO toolkit shared libraries\n" )
ELSE()
	MESSAGE( STATUS "Could not find ACADO toolkit shared libraries\n" )
ENDIF()

IF( VERBOSE )
	MESSAGE( STATUS "${ACADO_SHARED_LIBRARIES}\n" )
ENDIF()

#
# qpOASES embedded source files and include folders
# TODO: checks
#
SET( ACADO_QPOASES_EMBEDDED_SOURCES $ENV{ACADO_ENV_QPOASES_EMBEDDED_SOURCES} )
SET( ACADO_QPOASES_EMBEDDED_INC_DIRS $ENV{ACADO_ENV_QPOASES_EMBEDDED_INC_DIRS} )


#
# And finally set found flag...
#
IF( ACADO_INCLUDE_DIRS_FOUND AND ACADO_LIBRARY_DIRS_FOUND 
		AND ACADO_STATIC_LIBS_FOUND AND ACADO_SHARED_LIBS_FOUND )
	SET( ACADO_FOUND TRUE )
ENDIF()

MESSAGE( STATUS "********************************************************************************" )

