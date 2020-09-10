CC := gcc

# OPTIMIZED VS. DEBUGGING -- CHOOSE ONE
#CFLAGS := -I../ -Wall -ggdb -Wuninitialized 
#CFLAGS := -I../ -Wall -ggdb
#CFLAGS := -I../ -Wall -O1 
CFLAGS := -I../ -Wall  -ggdb

FF := /usr/bin/gfortran
#FF := /usr/local/gcc4/bin/gfortran
#FF := gfortran

# OPTIMIZED VS. DEBUGGING -- CHOOSE ONE
#FFLAGS := -fno-second-underscore -I -finit-local-zero -O1 -Wuninitialized
#FFLAGS := -fno-second-underscore -I -ggdb
#FFLAGS := -fno-second-underscore -I -ggdb -finit-local-zero -O1 -Wuninitialized
#FFLAGS := -fno-second-underscore -I -finit-local-zero -O1 -Wuninitialized -ggdb
#FFLAGS := -fno-second-underscore -ggdb
#FFLAGS := -fno-second-underscore -O1 
FFLAGS := -fno-second-underscore  -ggdb

CENTURY_LIB = libcentmc2.a

CENTURY_OBJS = adjlig.o agdrat.o anerob.o annacc.o atanf.o bgdrat.o \
	calciv.o candec.o cen_warmstart.o ckdata.o cmplig.o \
	co2eff.o crop.o cropin.o cropmix.o csa_detiv.o csa_main.o csched.o \
	cultiv.o cutrtn.o cycle.o daylen.o declig.o decomp.o \
	dedrem.o droot.o dshoot.o eachyr.o eq_test.o erosn.o esched.o \
	falstd.o faterr.o firein.o firrtn.o fixin.o \
	floclr.o flowUtilities.o frem.o frespr.o fsfunc.o getlen.o \
	gpdf.o grazin.o grem.o grochk.o growth.o h2olos.o harvst.o \
	init_climate.o init_lat.o init_soils.o inprac.o irrigt.o \
	killrt.o lacalc1.o lacalc2.o laprod.o leach.o line.o litburn.o \
	litdec.o livrem.o message.o mnracc.o nutrlm.o partit.o \
	passvars.o pevap.o potcrp.o potfor.o pprdwc.o prcgrw.o predec.o \
	prelim.o pschem.o ramp.o readblk.o repeat_event.o \
	respir.o restrp.o rtimp.o savarp.o scale.o scale4.o schedl.o \
	simsom.o sitein.o soilos.o somdec.o store_event.o sumcar.o \
	tcalc.o treein.o treemix.o trees.o tremin.o update_event.o \
	update_sched.o veg_change.o vetofix.o wdeath.o \
	weathr.o woodec.o wrtbin.o wthini.o

%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o : %.F
	$(FF) $(FFLAGS) -c $< -o $@


$(CENTURY_LIB): $(CENTURY_OBJS)
	rm -f $(CENTURY_LIB)
	ar q $(CENTURY_LIB) $(CENTURY_OBJS)
	ranlib $(CENTURY_LIB)

clean:
	rm -f $(CENTURY_OBJS) $(CENTURY_LIB)
