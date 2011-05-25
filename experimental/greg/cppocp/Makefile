PROJ = cppocp

Q ?= @

UNAME := $(shell uname)

CPP_SRC += \
	Ode.cpp \
	Ocp.cpp \
	OcpMultipleShooting.cpp \
	rocket.cpp
#	glider.cpp

LDFLAGS = \
	-rdynamic \
	-L/home/greg/casadi/build/lib \
	-L/usr/local/lib/coin \
	-L/usr/local/lib/coin/ThirdParty \
	-lipopt_interface \
	-lcasadi \
	-lipopt \
	-lpthread \
	-llapack \
	-lblas \
	-ldl \
	-lcoinhsl \
	-lcoinmetis \
	-Wl \
	-ldmumps_seq-4.9.2 \

INCLUDES = \
	-I/home/greg/casadi \
	-I/usr/local/include \
	-I/home/greg/casadi/external_packages/CSparse/Include

OPTFLAGS = -O3

OBJ = $(CPP_SRC:%.cpp=%.o) 

WARNINGFLAGS = #-Wall -Wextra -Wshadow -Werror

CXXFLAGS = $(WARNINGFLAGS) $(INCLUDES) $(OPTFLAGS) # -std=gnu++0x

CXXFLAGS += -fPIC -fopenmp -O3 -DNDEBUG

#CXX_DEFINES = -D_USE_MATH_DEFINES

CXX = g++

.PHONY: clean


$(PROJ): $(OBJ)
	@echo LD $@
	$(Q)$(CXX) $(CXXFLAGS) $(OBJ) $(LDFLAGS) -o $@

%.o : %.cpp
	@echo CXX $@
	$(Q)$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(PROJ)
	rm -f $(OBJ)
