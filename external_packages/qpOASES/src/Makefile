##
##	This file is part of qpOASES.
##
##	qpOASES -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2007-2012 by Hans Joachim Ferreau, Andreas Potschka,
##	Christian Kirches et al. All rights reserved.
##
##	qpOASES is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qpOASES is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qpOASES; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##



##
##	Filename:  src/Makefile
##	Author:    Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
##	Version:   3.0beta
##	Date:      2007-2012
##

SRCDIR = .

# select your operating system here!
include make_linux.mk
#include make_windows.mk


##
##	flags
##

IDIR        =  ../include
EXTRAS_DIR  =  extras

IFLAGS      =  -I. \
               -I${IDIR}


QPOASES_OBJECTS = \
	SQProblem.${OBJEXT} \
	QProblem.${OBJEXT} \
	QProblemB.${OBJEXT} \
	Bounds.${OBJEXT} \
	Constraints.${OBJEXT} \
	SubjectTo.${OBJEXT} \
	Indexlist.${OBJEXT} \
	Flipper.${OBJEXT} \
	Utils.${OBJEXT} \
	Options.${OBJEXT} \
	Matrices.${OBJEXT} \
	MessageHandling.${OBJEXT}


QPOASES_EXTRAS_OBJECTS = \
	SolutionAnalysis.${OBJEXT} \
	OQPinterface.${OBJEXT}



##
##	targets
##

all: libqpOASES.${LIBEXT} libqpOASESextras.${LIBEXT} ${LIB_LAPACK} ${LIB_BLAS}

withextras: libqpOASESextras.${LIBEXT}


libqpOASES.${LIBEXT}: ${QPOASES_OBJECTS}
	${AR} r $@ $?

libqpOASESextras.${LIBEXT}: ${QPOASES_OBJECTS} ${QPOASES_EXTRAS_OBJECTS}
	${AR} r $@ $?


SQProblem.${OBJEXT}: \
	SQProblem.cpp \
	${IDIR}/qpOASES/SQProblem.hpp \
	${IDIR}/qpOASES/QProblem.hpp \
	${IDIR}/qpOASES/QProblemB.hpp \
	${IDIR}/qpOASES/Bounds.hpp \
	${IDIR}/qpOASES/Constraints.hpp \
	${IDIR}/qpOASES/SubjectTo.hpp \
	${IDIR}/qpOASES/Indexlist.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/Matrices.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} SQProblem.cpp

QProblem.${OBJEXT}: \
	QProblem.cpp \
	${IDIR}/qpOASES/QProblem.hpp \
	${IDIR}/qpOASES/QProblemB.hpp \
	${IDIR}/qpOASES/Bounds.hpp \
	${IDIR}/qpOASES/Constraints.hpp \
	${IDIR}/qpOASES/SubjectTo.hpp \
	${IDIR}/qpOASES/Indexlist.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/Matrices.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} QProblem.cpp

QProblemB.${OBJEXT}: \
	QProblemB.cpp \
	${IDIR}/qpOASES/QProblemB.hpp \
	${IDIR}/qpOASES/Bounds.hpp \
	${IDIR}/qpOASES/Constraints.hpp \
	${IDIR}/qpOASES/SubjectTo.hpp \
	${IDIR}/qpOASES/Indexlist.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/Matrices.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} QProblemB.cpp

Bounds.${OBJEXT}: Bounds.cpp \
	${IDIR}/qpOASES/Bounds.hpp \
	${IDIR}/qpOASES/SubjectTo.hpp \
	${IDIR}/qpOASES/Indexlist.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} Bounds.cpp

Constraints.${OBJEXT}: \
	Constraints.cpp \
	${IDIR}/qpOASES/Constraints.hpp \
	${IDIR}/qpOASES/SubjectTo.hpp \
	${IDIR}/qpOASES/Indexlist.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} Constraints.cpp

SubjectTo.${OBJEXT}: \
	SubjectTo.cpp \
	${IDIR}/qpOASES/SubjectTo.hpp \
	${IDIR}/qpOASES/Indexlist.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} SubjectTo.cpp

Indexlist.${OBJEXT}: \
	Indexlist.cpp \
	${IDIR}/qpOASES/Indexlist.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} Indexlist.cpp

Flipper.${OBJEXT}: \
	Flipper.cpp \
	${IDIR}/qpOASES/Flipper.hpp \
	${IDIR}/qpOASES/Bounds.hpp \
	${IDIR}/qpOASES/Constraints.hpp \
	${IDIR}/qpOASES/SubjectTo.hpp \
	${IDIR}/qpOASES/Indexlist.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} Flipper.cpp

Utils.${OBJEXT}: \
	Utils.cpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} Utils.cpp

Options.${OBJEXT}: \
	Options.cpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} Options.cpp

Matrices.${OBJEXT}: \
	Matrices.cpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} Matrices.cpp

BLASReplacement.${OBJEXT}: \
	BLASReplacement.cpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} BLASReplacement.cpp

LAPACKReplacement.${OBJEXT}: \
	LAPACKReplacement.cpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} LAPACKReplacement.cpp

MessageHandling.${OBJEXT}: \
	MessageHandling.cpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} MessageHandling.cpp

SolutionAnalysis.${OBJEXT}: \
	${EXTRAS_DIR}/SolutionAnalysis.cpp \
	${IDIR}/qpOASES/${EXTRAS_DIR}/SolutionAnalysis.hpp \
	${IDIR}/qpOASES/SQProblem.hpp \
	${IDIR}/qpOASES/QProblem.hpp \
	${IDIR}/qpOASES/QProblemB.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/Matrices.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} ${EXTRAS_DIR}/SolutionAnalysis.cpp

OQPinterface.${OBJEXT}: \
	${EXTRAS_DIR}/OQPinterface.cpp \
	${IDIR}/qpOASES/${EXTRAS_DIR}/OQPinterface.hpp \
	${IDIR}/qpOASES/QProblem.hpp \
	${IDIR}/qpOASES/QProblemB.hpp \
	${IDIR}/qpOASES/Utils.hpp \
	${IDIR}/qpOASES/Constants.hpp \
	${IDIR}/qpOASES/Types.hpp \
	${IDIR}/qpOASES/Options.hpp \
	${IDIR}/qpOASES/Matrices.hpp \
	${IDIR}/qpOASES/MessageHandling.hpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} ${EXTRAS_DIR}/OQPinterface.cpp

clean:
	${RM} -f *.${OBJEXT} *.${LIBEXT}

clobber: clean


%.${OBJEXT}: %.cpp
	@echo "Creating" $@
	${CPP} ${DEF_TARGET} ${IFLAGS} ${CPPFLAGS} -c $<

shared:	${QPOASES_OBJECTS}
	${CPP} -shared -Wl,-soname,libqpOASES.so -o libqpOASES.so ${QPOASES_OBJECTS} -lm -llapack -lblas

##
##	end of file
##
