
PRJNAME := saruman

# Source files here
EXECUTABLE	:= $(PRJNAME)

# Cuda source files (compiled with cudacc)
CUFILES		:= cuda-functions.cu

# C/C++ source files (compiled with gcc / c++)
CFILES		:= saruman.c util.c hash.c

# Cuda AND C/C++ build-dependencies
CU_AND_C_DEPS := cuda-functions.h

# Cuda build-dependencies
CU_DEPS     := $(CU_AND_C_DEPS)

# C/C++ build-dependencies
C_DEPS      := saruman.c saruman.h util.c util.h hash.c hash.h $(CU_AND_C_DEPS)

################################################################################
# Rules and targets
################################################################################

.SUFFIXES : .cu .cu_dbg_o .c_dbg_o .cpp_dbg_o .cu_rel_o .c_rel_o .cpp_rel_o .cubin

# Basic directory setup for Project
ROOTDIR    := .
SRCDIR     := $(ROOTDIR)
ROOTBINDIR := ../build/bin
BINDIR     := $(ROOTBINDIR)/linux
ROOTOBJDIR := ../build/obj
CUDAVERSION :=

# Basic directory setup for CUDA
CUDA_INSTALL_PATH := /mnt/misc/cuda/10.1/
CUDA_SDK_DIR := /mnt/misc/cuda/10.1/
LIBDIR     := $(CUDA_SDK_DIR)/lib64
COMMONDIR  := $(CUDA_SDK_DIR)/samples/common

# Switching device-emulation-mode (emu := 1 is switched on)
emu := 0
dbg := 0

# Switching verbose-mode (verbose := 1 is switched on for more compile-information)
verbose := 1

# architecture/version-flags. Here for cc1.1
#SMVERSIONFLAGS	:= -arch sm_20 #-gencode arch=compute_20,code=sm_20
	
# Includes
INCLUDES  += -I. -I$(CUDA_INSTALL_PATH)/samples/common/inc/ -I$(CUDA_INSTALL_PATH)/include -I$(COMMONDIR)/inc -I /usr/local/lib -I /homes/tjakobi/local/include/ -I /homes/tjakobi/local/lib/

# Compilers
NVCC       := /mnt/misc/cuda/10.1/bin/nvcc
CXX        := g++ 
CC         := gcc 
LINK       := gcc -fPIC $(INCLUDES)

# Warning flags
CXXWARN_FLAGS := \
	-W -Wall \
	-Wimplicit \
	-Wswitch \
	-Wformat \
	-Wchar-subscripts \
	-Wparentheses \
	-Wmultichar \
	-Wtrigraphs \
	-Wpointer-arith \
	-Wcast-align \
	-Wreturn-type \
	-Wno-unused-function \
	$(SPACE)

CWARN_FLAGS := $(CXXWARN_FLAGS) \
	-Wstrict-prototypes \
	-Wmissing-prototypes \
	-Wmissing-declarations \
	-Wnested-externs \
	-Wmain \

# Compiler-specific flags
NVCCFLAGS := --ptxas-options -v #--maxrregcount=16
CXXFLAGS  := $(CXXWARN_FLAGS)
CFLAGS    := $(CWARN_FLAGS) -std=c99 -pedantic # -DHASH_EMIT_KEYS=3

# Common flags
COMMONFLAGS = $(INCLUDES)

# Debug/release configuration
ifeq ($(dbg),1)
	COMMONFLAGS += -g  -DHASH_FUNCTION=HASH_FNV
	NVCCFLAGS   += -D_DEBUG -G
	BINSUBDIR   := debug
	LIBSUFFIX   := x86_64
else 
	COMMONFLAGS += -DHASH_FUNCTION=HASH_FNV 
	BINSUBDIR   := release
	LIBSUFFIX   := _x86_64
	NVCC_FLAGS  += --compiler-options -fno-strict-aliasing
	CXXFLAGS  	+= -fno-strict-aliasing 
	CFLAGS  	+= -fno-strict-aliasing -O3 -pipe -msse3 
endif

# append optional arch/SM version flags (such as -arch sm_11)
NVCCFLAGS += $(SMVERSIONFLAGS)

# architecture flag for cubin build
CUBIN_ARCH_FLAG := -m64

# detect if 32 bit or 64 bit system
HP_64 =	$(shell uname -m | grep 64)

# OpenGL is used or not (if it is used, then it is necessary to include GLEW)
OPENGLLIB := -lGL -lGLU -lglut
ifeq ($(USEGLLIB),1)
	ifeq "$(strip $(HP_64))" ""
		OPENGLLIB += -lGLEW
	else
		OPENGLLIB += -lGLEW_x86_64
	endif

	CUBIN_ARCH_FLAG := -m64
endif

ifeq ($(USEPARAMGL),1)
	PARAMGLLIB := -lparamgl$(LIBSUFFIX)
endif

ifeq ($(USECUDPP), 1)
	ifeq "$(strip $(HP_64))" ""
		CUDPPLIB := -lcudpp
	else
		CUDPPLIB := -lcudpp64
	endif

	CUDPPLIB := $(CUDPPLIB)$(LIBSUFFIX)

	ifeq ($(emu), 1)
		CUDPPLIB := $(CUDPPLIB)_emu
	endif
endif

# Libs

LIB       := -L/mnt/misc/cuda/10.1/targets/x86_64-linux/lib/stubs -L$(CUDA_INSTALL_PATH)/lib64 -L$(LIBDIR) -L$(CUDA_SDK_DIR)/C/lib -L$(COMMONDIR)/lib -L /homes/tjakobi/local/lib/ -Wl,-Bstatic -largtable2 -Wl,-Bdynamic -lpthread -lcuda -lm 
 
ifeq ($(USEDRVAPI),1)
   LIB += -lcuda -lcudart ${OPENGLLIB} $(PARAMGLLIB) $(CUDPPLIB) ${LIB} 
else
   LIB += -lcudart $(PARAMGLLIB) $(CUDPPLIB) ${LIB}
endif

ifeq ($(USECUFFT),1)
  ifeq ($(emu),1)
    LIB += -lcufftemu
  else
    LIB += -lcufft
  endif
endif

ifeq ($(USECUBLAS),1)
  ifeq ($(emu),1)
    LIB += -lcublasemu
  else
    LIB += -lcublas
  endif
endif

# Lib/exe configuration
ifneq ($(STATIC_LIB),)
	TARGETDIR := $(LIBDIR)
	TARGET   := $(subst .a,$(LIBSUFFIX).a,$(LIBDIR)/$(STATIC_LIB))
	LINKLINE  = ar qv $(TARGET) $(OBJS) 
else
	# LIB += -lcutil$(LIBSUFFIX)
	# Device emulation configuration
	ifeq ($(emu), 1)
		NVCCFLAGS   += -deviceemu
		CUDACCFLAGS += 
		BINSUBDIR   := emu$(BINSUBDIR)
		# consistency, makes developing easier
		CXXFLAGS		+= -D__DEVICE_EMULATION__
		CFLAGS			+= -D__DEVICE_EMULATION__
	endif
	TARGETDIR := $(BINDIR)/$(BINSUBDIR)
	TARGET    := $(TARGETDIR)/$(EXECUTABLE)
	LINKLINE  = $(LINK) -o $(TARGET) $(OBJS) $(LIB)
endif

# check if verbose 
ifeq ($(verbose), 1)
	VERBOSE :=
else
	VERBOSE := @
endif

################################################################################
# Check for input flags and set compiler flags appropriately
################################################################################
ifeq ($(fastmath), 1)
	NVCCFLAGS += -use_fast_math
endif

ifeq ($(keep), 1)
	NVCCFLAGS += -keep
	NVCC_KEEP_CLEAN := *.i* *.cubin *.cu.c *.cudafe* *.fatbin.c *.ptx
endif

ifdef maxregisters
	NVCCFLAGS += -maxrregcount $(maxregisters)
endif

# Add cudacc flags
NVCCFLAGS += $(CUDACCFLAGS)

# Add common flags
NVCCFLAGS += $(COMMONFLAGS)
CXXFLAGS  += $(COMMONFLAGS)
CFLAGS    += $(COMMONFLAGS)

ifeq ($(nvcc_warn_verbose),1)
	NVCCFLAGS += $(addprefix --compiler-options ,$(CXXWARN_FLAGS)) 
	NVCCFLAGS += --compiler-options -fno-strict-aliasing
endif

################################################################################
# Set up object files
################################################################################
OBJDIR := $(ROOTOBJDIR)/$(BINSUBDIR)
OBJS :=  $(patsubst %.cpp,$(OBJDIR)/%.cpp_o,$(notdir $(CCFILES)))
OBJS +=  $(patsubst %.c,$(OBJDIR)/%.c_o,$(notdir $(CFILES)))
OBJS +=  $(patsubst %.cu,$(OBJDIR)/%.cu_o,$(notdir $(CUFILES)))

################################################################################
# Set up cubin files
################################################################################
CUBINDIR := $(SRCDIR)/data
CUBINS +=  $(patsubst %.cu,$(CUBINDIR)/%.cubin,$(notdir $(CUBINFILES)))

################################################################################
# Rules
################################################################################

run : 
	$(TARGET)

cleanrun : clean all

makeAllAndRun : all
	@echo "running"
	$(TARGET)

all : $(TARGET) Makefile $(OBJS)
	@echo "make $@"
	@echo "...make all done."

$(TARGET): makedirectories $(OBJS) $(CUBINS) Makefile
	@echo "make $@ (target)"
	$(VERBOSE)$(LINKLINE)

$(OBJDIR)/%.c_o : $(SRCDIR)/%.c $(C_DEPS)
	@echo "make $@"
	$(VERBOSE)$(CC) $(CFLAGS) -o $@ -c $<

$(OBJDIR)/%.cpp_o : $(SRCDIR)/%.cpp $(C_DEPS)
	@echo "make $@"
	$(VERBOSE)$(CXX) $(CXXFLAGS) -o $@ -c $<

$(OBJDIR)/%.cu_o : $(SRCDIR)/%.cu $(CU_DEPS)
	@echo "make $@"
	$(VERBOSE)$(NVCC) -o $@ -c $< $(NVCCFLAGS)

$(CUBINDIR)/%.cubin : cubindirectory $(SRCDIR)/%.cu
	@echo "make $@"
	$(VERBOSE)$(NVCC) -m32 -o $@ -cubin $< $(NVCCFLAGS)

cubindirectory:
	@echo "make $@"
	$(VERBOSE)mkdir -p $(CUBINDIR)

makedirectories:
	@echo "make $@"
	$(VERBOSE)mkdir -p $(OBJDIR)
	$(VERBOSE)mkdir -p $(TARGETDIR)

.PHONY : clean
clean :
	@echo "make clean"
	$(VERBOSE)rm -f $(OBJS)
	$(VERBOSE)rm -f $(CUBINS)
	$(VERBOSE)rm -f $(TARGET)
	@echo "...make clean done."

.PHONY : cleanroot
cleanall :
	@echo "make cleanall"
	$(VERBOSE)rm -rdf $(ROOTBINDIR)
	$(VERBOSE)rm -rdf $(ROOTOBJDIR)	
	@echo "...make cleanall done."
