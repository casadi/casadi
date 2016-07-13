"""
 blockSQP -- Sequential quadratic programming for problems with
             block-diagonal Hessian matrix.
 Copyright (C) 2012-2015 by Dennis Janka <dennis.janka@iwr.uni-heidelberg.de>

 Licensed under the zlib license. See LICENSE for more details.
"""

#
# Create a single header file from
#
# 	blocksqp_matrix.hpp
# 	blocksqp_problemspec.hpp
#
# Then modify
#
# 	blocksqp_matrix.cpp
# 	blocksqp_problemspec.cpp
#
# accordingly and create a makefile that compile these into a library
#
# 	libblocksqp_min.so
#
# Make a tar-file of the whole package and copy to vplan/dist
import os

incDir = "include/"
srcDir = "src/"
minDir = "blockSQP_min/"
if not os.path.exists(minDir):
    os.makedirs(minDir)

#
# 1.) CREATE NEW HEADER FILE
#
# read blocksqp_defs.hpp
myFile = open( incDir + "blocksqp_defs.hpp", "r" )
defsLines = myFile.readlines()
myFile.close()

# read blocksqp_matrix.hpp
myFile = open( incDir + "blocksqp_matrix.hpp", "r" )
matrixLines = myFile.readlines()
myFile.close()

# read blocksqp_problemspec.hpp
myFile = open( incDir + "blocksqp_problemspec.hpp", "r" )
probLines = myFile.readlines()
myFile.close()

# insert blocksqp_defs.hpp where it is included in blocksqp_matrix.hpp
k = 0
for line in matrixLines:
    if "#include" in line and "blocksqp_defs" in line:
        defsLines.reverse()
        for insertLine in defsLines:
            matrixLines.insert( k, insertLine )
        matrixLines.remove( line )
        break
    k += 1

# insert blocksqp_matrix.hpp where it is included in blocksqp_problemspec.hpp
k = 0
for line in probLines:
    if "#include" in line and "blocksqp_matrix" in line:
        matrixLines.reverse()
        for insertLine in matrixLines:
            probLines.insert( k, insertLine )
        probLines.remove( line )
        break
    k += 1

# write combined header files to new file
myFile = open( minDir + "blocksqp_min.hpp", "w" )
for line in probLines:
    myFile.writelines( line )
myFile.close()


#
# 2.) CREATE NEW CPP FILES
#
# read blocksqp_matrix.cpp and replace include
myFile = open( srcDir + "blocksqp_matrix.cpp", "r" )
lines = myFile.readlines()
myFile.close()
k=0
for line in lines:
    if "#include" in line and "blocksqp_matrix" in line:
        lines.insert( k, "#include \"blocksqp_min.hpp\"\n" )
        lines.remove( line )
        break
    k += 1

# write new cpp file
myFile = open( minDir + "blocksqp_matrix.cpp", "w" )
for line in lines:
    myFile.writelines( line )
myFile.close()

# read blocksqp_problemspec.cpp and replace include
myFile = open( srcDir + "blocksqp_problemspec.cpp", "r" )
lines = myFile.readlines()
myFile.close()
k=0
for line in lines:
    if "#include" in line and "blocksqp_problemspec" in line:
        lines.insert( k, "#include \"blocksqp_min.hpp\"\n" )
        lines.remove( line )
        break
    k += 1

# write new cpp file
myFile = open( minDir + "blocksqp_problemspec.cpp", "w" )
for line in lines:
    myFile.writelines( line )
myFile.close()


#
# 3.) CREATE MAKEFILE
#
myFile = open( minDir + "makefile", "w" )
myFile.writelines( "# Makefile for minimum version\n" )
myFile.writelines( "OBJECTS = blocksqp_matrix.o \\\n" )
myFile.writelines( "\t\tblocksqp_problemspec.o\n\n" )
myFile.writelines( "OPTIONS = -g -O0 -fPIC -Wno-deprecated -Wno-write-strings\n\n" )
myFile.writelines( "library: $(OBJECTS)\n" )
myFile.writelines( "\tg++ -shared -o libblockSQP_min.so $(OBJECTS)\n\n" )
myFile.writelines( "%.o: %.cpp\n" )
myFile.writelines( "\tg++ -c $(OPTIONS) -o $@ $<\n\n" )
myFile.writelines( "clean:\n" )
myFile.writelines( "\trm -rf *.o libblockSQP_min.so\n" )
myFile.close()


#
# 4.) CREATE TAR FILE AND COPY TO vplan/dist
#
os.system( "tar -czf  blockSQP_min.tar.gz blockSQP_min/" )
os.system( "mv blockSQP_min.tar.gz $VPLANROOT/dist" )
os.system( "rm -rf blockSQP_min" )












