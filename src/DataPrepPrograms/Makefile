CPP := /usr/bin/g++
CPPFLAGS := -I/opt/local/include -Wall -ggdb -DUNIX_MAPSS -O1 -Wuninitialized
LD := /usr/bin/gfortran
LDFLAGS := -fno-second-underscore -lstdc++

NETCDF_LIB := -L/opt/local/lib/ -lnetcdf

MAKENCFILE := makeNCfile

MAKENCFILE_OBJS := makeNCfile.o

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

$(MAKENCFILE): $(MAKENCFILE_OBJS)  
	$(LD) $(LDFLAGS) -o makeNCfile $(MAKENCFILE_OBJS) $(NETCDF_LIB)

MERGENCDATA := mergeNCdata

MERGENCDATA_OBJS := mergeNCdata.o

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

$(MERGENCDATA): $(MERGENCDATA_OBJS)
	$(LD) $(LDFLAGS) -o mergeNCdata $(MERGENCDATA_OBJS) $(NETCDF_LIB)

MAKEEQ := makeEQ

MAKEEQ_OBJS := makeEQ.o

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

$(MAKEEQ): $(MAKEEQ_OBJS)
	$(LD) $(LDFLAGS) -o makeEQ $(MAKEEQ_OBJS) $(NETCDF_LIB)

MAKESPINUP := makeSpinup

MAKESPINUP_OBJS := makeSpinup.o

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

$(MAKESPINUP): $(MAKESPINUP_OBJS)
	$(LD) $(LDFLAGS) -o makeSpinup $(MAKESPINUP_OBJS) $(NETCDF_LIB)

ADD_ATTRIBUTE := add_attribute

ADD_ATTRIBUTE_OBJS := add_attribute.o

%.o : %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@

$(ADD_ATTRIBUTE): $(ADD_ATTRIBUTE_OBJS)
	$(LD) $(LDFLAGS) -o add_attribute $(ADD_ATTRIBUTE_OBJS) $(NETCDF_LIB)

clean:
	rm -f $(MAKENCFILE_OBJS) $(MAKENCFILE)
	rm -f $(MERGENCDATA_OBJS) $(MERGENCDATA)
	rm -f $(MAKEEQ_OBJS) $(MAKEEQ)
	rm -f $(MAKESPINUP_OBJS) $(MAKESPINUP)
	rm -f $(ADD_ATTRIBUTE_OBJS) $(ADD_ATTRIBUTE)

