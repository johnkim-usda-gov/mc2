#CPP := /usr/bin/g++
CPP := /usr/bin/g++

# PRODUCTION v. DEBUGGING
#CPPFLAGS := -I/opt/local/include -Wall -DUNIX_MAPSS -Wuninitialized -O
#CPPFLAGS := -I/opt/local/include -Wall -ggdb -DUNIX_MAPSS 
#CPPFLAGS := -DUNIX_MAPSS -O1
CPPFLAGS := -DUNIX_MAPSS -ggdb

LD := /usr/bin/gfortran
#LD := /usr/local/gcc4/bin/gfortran
#LD := gfortran

LDFLAGS := -fno-second-underscore -lstdc++

#NETCDF_LIB := -L/opt/local/lib/ -lnetcdf
NETCDF_LIB := -L/usr/local/lib/ -lnetcdf

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
  commandFile.o

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
