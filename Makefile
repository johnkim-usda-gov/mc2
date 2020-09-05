# CPP := /usr/bin/g++ gcc version 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.9.00)
# CPP := /usr/local/gfortran/bin/g++ gcc version 4.8.2 (GCC)
CPP := /usr/bin/g++
CPPFLAGS := -I/opt/local/include -Wall -ggdb -DUNIX_MAPSS -O1 -Wuninitialized
LD := /usr/bin/gfortran
LDFLAGS := -fno-second-underscore -lstdc++

NETCDF_LIB := -L/opt/local/lib/ -lnetcdf

CENTURY_DIR := ./century/
CENTURY_LIB := $(CENTURY_DIR)/libcentmc2.a

MC2 := mc2

MC2_OBJS := \
  MC2.o \
  ProcessModel.o \
  ScienceFcns.o \
  MAPSSbiogeographyModel.o \
  MAPSSfcns.o \
  CENTURY.o \
  MCfire.o \
  MCbiogeog.o \
  mc2_secondary.o \
  commandFile.o \
  PNVmodel.o

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

$(MC2): $(MC2_OBJS)  
	$(LD) $(LDFLAGS) -o mc2 $(MC2_OBJS) $(CENTURY_LIB) $(NETCDF_LIB)

.PHONY: dummy clean realclean


clean:
	rm -f $(MC2_OBJS) $(MC2)

realclean:
	make clean
	(cd $(CENTURY_DIR); make clean)
