##
## blockSQP -- Sequential quadratic programming for problems with
##             block-diagonal Hessian matrix.
## Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>
##
## Licensed under the zlib license. See LICENSE for more details.
##

########################################################################
#                User configuration: qpOASES settings                  #
########################################################################

# Location of qpOASES header files and shared library
QPOASESDIR = ~/path/to/qpOASES
QPOASESINCLUDE = $(QPOASESDIR)/include
QPOASESLIBDIR = $(QPOASESDIR)/bin

########################################################################
#                      End of user configuration                       #
########################################################################

INCLUDEDIR = include
SRCDIR = src
EXDIR = examples
LIBDIR = lib
OBJDIR = $(LIBDIR)/obj

OPTIONS = -g -O0 -fPIC -I $(INCLUDEDIR) -I $(QPOASESINCLUDE) -Wno-deprecated -Wno-write-strings -Wall
#OPTIONS = -g -O3 -fPIC -I $(INCLUDEDIR) -I $(QPOASESINCLUDE) -Wno-deprecated -Wno-write-strings -Wall

EXAMPLES = $(EXDIR)/example1

OBJECTS = $(OBJDIR)/blocksqp_matrix.o \
		$(OBJDIR)/blocksqp_problemspec.o \
		$(OBJDIR)/blocksqp_general_purpose.o \
		$(OBJDIR)/blocksqp_glob.o \
		$(OBJDIR)/blocksqp_hess.o \
		$(OBJDIR)/blocksqp_iter.o \
		$(OBJDIR)/blocksqp_main.o \
		$(OBJDIR)/blocksqp_options.o \
		$(OBJDIR)/blocksqp_qp.o \
		$(OBJDIR)/blocksqp_restoration.o \
		$(OBJDIR)/blocksqp_stats.o

LIBS       = -L $(CURDIR)/$(LIBDIR) -Xlinker -rpath -Xlinker $(CURDIR)/$(LIBDIR) \
             -L $(QPOASESLIBDIR) -Xlinker -rpath -Xlinker $(QPOASESLIBDIR) \
             -L /usr/local/lib -Xlinker -rpath -Xlinker /usr/local/lib \
             -lblockSQP \
             -lqpOASES \
             -llapack

all: library examples

library: $(OBJECTS) | $(LIBDIR)
	g++ -shared -o $(LIBDIR)/libblockSQP.so $(OBJECTS)

min: $(OBJECTS) | $(LIBDIR)
	g++ -shared -o $(LIBDIR)/libblockSQP_min.so $(OBJECTS)

examples: $(EXAMPLES)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp makefile | $(OBJDIR)
	g++ -c $(OPTIONS) -o $@ $<

$(EXDIR)/%: $(EXDIR)/%.cpp makefile
	g++ $(OPTIONS) -o $@ $< $(LIBS)

$(OBJDIR): | $(LIBDIR)
	mkdir $(OBJDIR)

$(LIBDIR):
	mkdir $(LIBDIR)

.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(LIBDIR) $(EXAMPLES)
