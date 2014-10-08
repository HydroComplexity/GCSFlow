# Make file for GCS-Flow model
# Copyright 2014 Phong Le <levuvietphong@gmail.com>
# All Rights Reserved.

# Executable name
TARGET = gcsflow

# Project Name
BINDIR = bin

# Source, Include, Object, and Library Directories
SRCDIR = src
INCDIR = include
OBJDIR = src/obj
LIBDIR = lib

# Compiler, Flags, and Library
NVCC = nvcc
NVCCFLAGS = -O3 -arch=sm_20 -I/$(INCDIR)
GXX = g++
GXXFLAGS = -O3 -Wall
LIBS = -lnetcdf
CUSOURCES = $(shell cd $(SRCDIR) && ls *.cu)
CCSOURCES = $(shell cd $(SRCDIR) && ls *.cc)

# CUDA source code and dependency
CUDEPS = $(patsubst %,$(INCDIR)/%,$(CUSOURCES:.cu=.h))
CUOBJ = $(patsubst %,$(OBJDIR)/%,$(CUSOURCES:.cu=.cu.o))

# CPP source code and dependency
CCDEPS = $(patsubst %,$(INCDIR)/%,$(CCSOURCES:.cc=.h))
CCOBJ = $(patsubst %,$(OBJDIR)/%,$(CCSOURCES:.cc=.cc.o))

all: build print $(TARGET)

$(OBJDIR)/%.cu.o: $(SRCDIR)/%.cu $(CUDEPS)
	$(NVCC) $(NVCCFLAGS) -dc -o $@ $<

$(OBJDIR)/%.cc.o: $(SRCDIR)/%.cc $(CCDEPS)
	$(GXX) $(GXXFLAGS) -c -o $@ $<

lib_cuda.a: $(CUOBJ)
	$(NVCC) $(NVCCFLAGS) -lib $^ -o $(BINDIR)/$@

$(TARGET): lib_cuda.a $(CCOBJ)
	@echo ""
	@echo "Linking..."
	$(NVCC) $(NVCCFLAGS) $(OBJDIR)/*cc.o $(BINDIR)/$< -o $(BINDIR)/$@ $(LIBS)
	@echo ""
	@echo "Compilation file "$@" successfully. DONE!!!!"
	@echo ""

# Clean-up executables and object files
.PHONY: clean

print:
	@echo ""
	@echo "Compling objects..."

# Check for the existence of sub directories.
# If not exist, then create these directories to dump data to.
# If exist do nothing
build:
	@echo ""
	@echo "Creating Directories..."
	if [ ! -d "$(LIBDIR)" ]; then mkdir $(LIBDIR);  fi
	if [ ! -d "$(OBJDIR)" ]; then mkdir $(OBJDIR);  fi
	if [ ! -d "$(BINDIR)" ]; then mkdir $(BINDIR);  fi

clean:
	rm -f $(BINDIR)/lib_cuda.a $(OBJDIR)/*.o *~ core $(INCDIR)/*~
