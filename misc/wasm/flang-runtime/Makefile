ROOT = $(abspath .)
BUILD = $(ROOT)/build
SOURCE = $(ROOT)/llvm-project
PREFIX = $(ROOT)

HOST = $(PREFIX)/host
WASM = $(PREFIX)/wasm

FLANG_BIN = $(BUILD)/bin/flang
FLANG_INCLUDE = $(BUILD)/include/flang

RUNTIME_SOURCES := $(wildcard $(SOURCE)/flang-rt/lib/runtime/*.cpp)
RUNTIME_SOURCES += $(SOURCE)/flang/lib/Decimal/decimal-to-binary.cpp
RUNTIME_SOURCES += $(SOURCE)/flang/lib/Decimal/binary-to-decimal.cpp
RUNTIME_OBJECTS = $(patsubst $(SOURCE)/%,$(BUILD)/%,$(RUNTIME_SOURCES:.cpp=.o))
RUNTIME_LIB = $(BUILD)/flang/runtime/libFortranRuntime.a
RUNTIME_CONFIG = $(BUILD)/tools/flang/runtime/config.h

FLANG_WASM_CMAKE_VARS := $(FLANG_WASM_CMAKE_VARS)

.PHONY: all
all: flang wasm-runtime

.PHONY: download
download: $(SOURCE)

.PHONY: wasm-runtime
wasm-runtime: $(RUNTIME_LIB)

$(SOURCE):
	git clone --single-branch -b wasm --depth=1 https://github.com/r-wasm/llvm-project

.PHONY: flang
flang: $(FLANG_BIN)
$(FLANG_BIN): $(SOURCE)
	@mkdir -p $(BUILD)
	cmake -G Ninja -S $(SOURCE)/llvm -B $(BUILD) \
	  -DCMAKE_INSTALL_PREFIX=$(HOST) \
	  -DCMAKE_BUILD_TYPE=MinSizeRel \
	  -DLLVM_DEFAULT_TARGET_TRIPLE="wasm32-unknown-emscripten" \
	  -DLLVM_TARGETS_TO_BUILD="WebAssembly" \
	  -DLLVM_ENABLE_PROJECTS="clang;flang;mlir" \
	  -DFLANG_INCLUDE_TESTS=OFF \
	  -DFLANG_ENABLE_FLANG_RT=OFF \
	  $(FLANG_WASM_CMAKE_VARS)
	TERM=dumb cmake --build $(BUILD)
	$(MAKE) wasm-runtime

RUNTIME_CXXFLAGS := $(RUNTIME_CXXFLAGS)
RUNTIME_CXXFLAGS += -I$(BUILD)/include -I$(BUILD)/tools/flang/runtime
RUNTIME_CXXFLAGS += -I$(SOURCE)/flang/include -I$(SOURCE)/flang-rt/include
RUNTIME_CXXFLAGS += -I$(SOURCE)/llvm/include
RUNTIME_CXXFLAGS += -DFLANG_LITTLE_ENDIAN
RUNTIME_CXXFLAGS += -fPIC -Wno-c++11-narrowing -fvisibility=hidden
RUNTIME_CXXFLAGS += -DFE_UNDERFLOW=0 -DFE_OVERFLOW=0 -DFE_INEXACT=0
RUNTIME_CXXFLAGS += -DFE_INVALID=0 -DFE_DIVBYZERO=0 -DFE_ALL_EXCEPT=0

$(RUNTIME_CONFIG):
	@mkdir -p $(@D)
	@printf '#ifndef FORTRAN_RUNTIME_CONFIG_H\n#define FORTRAN_RUNTIME_CONFIG_H\n#define HAVE_STRERROR_R 1\n#define HAVE_DECL_STRERROR_S 0\n/* #undef HAVE_BACKTRACE */\n#endif\n' > $@

$(RUNTIME_LIB): $(RUNTIME_OBJECTS)
	@mkdir -p $(@D)
	@rm -f $@
	emar -rcs $@ $^

$(BUILD)%.o : $(SOURCE)%.cpp $(RUNTIME_CONFIG)
	@mkdir -p $(@D)
	em++ $(RUNTIME_CXXFLAGS) -o $@ -c $<

.PHONY: install
install: $(FLANG_BIN) $(RUNTIME_LIB)
	mkdir -p $(HOST)/bin
	mkdir -p $(HOST)/include/flang
	mkdir -p $(WASM)/lib
	install -m 755 $(FLANG_BIN) $(HOST)/bin
	install -m 644 $(RUNTIME_LIB) $(WASM)/lib
	find $(FLANG_INCLUDE) -type f -exec install -m 644 {} $(HOST)/include/flang ';'

.PHONY: check
check:
	cmake --build $(BUILD) --target check-all

.PHONY: clean
clean:
	cmake --build $(BUILD) --target clean

.PHONY: clean-all
clean-all:
	rm -rf $(SOURCE) $(BUILD)
