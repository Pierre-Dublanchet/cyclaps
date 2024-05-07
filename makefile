FC=mpif90
FCFLAGS=-O5
NETCDF=/opt/homebrew/Cellar/netcdf/4.9.2_1
NETCDFF=/opt/homebrew/Cellar/netcdf-fortran/4.6.1
INC_SF=./src/special_functions/
GEOM=$(word 2, $(subst _, ,$(exec)))
KERNEL=$(word 3, $(subst _, ,$(exec)))
SLAW=$(word 4, $(subst _, ,$(exec)))
PP=$(word 5, $(subst _, ,$(exec)))
SC=$(word 6, $(subst _, ,$(exec)))

SRC_DIRS = $(sort $(dir ./src/ ./src/cfgio/ ./src/fft/ ./src/special_functions/))
vpath %.f90 $(SRC_DIRS)
vpath %.f $(SRC_DIRS)

OBJS_IN=string_conv_mod.o cfgio_mod.o

ifeq "$(KERNEL)" "freesurface"
	OBJS := mod_cst_fault_rns.o libcfgio.a
else
	OBJS := mod_cst_fault_rns.o libcfgio.a zfft1d.o factor.o pfactor.o fft235.o kernel.o
endif

ifeq "$(KERNEL)" "cr"
	OBJS :=  $(OBJS) libspecial_functions.a
endif

ifeq "$(PP)" "pnl"
	OBJS :=  $(OBJS) mod_cst_fault_rns_p.o
endif

ifeq "$(PP)" "press"
	OBJS :=  $(OBJS) mod_cst_fault_rns_p.o
endif

ifeq "$(SC)" "rs"
	OBJS :=  $(OBJS) mod_cst_fault_rns_p_rs.o
endif

ifneq "$(KERNEL)" "freesurface"
	ifeq "$(GEOM)" "2d"
		OBJS := $(OBJS) pzfft1d.o
	else
		OBJS := $(OBJS) pzfft2d.o
	endif
endif

OBJS := $(OBJS) $(exec).o
#OBJS_FILES = $(patsubst %.o,./obj/%.o,$(OBJS))
#OBJS_LIB_FILES = $(patsubst %.a,./lib/%.a,$(OBJS_FILES))

do: $(exec)
#	rm -f ./*.o ./*.a ./*.mod
	@echo $(NETCDF)
#@echo $(SRC_DIRS)
#@echo $(GEOM) 
#@echo $(KERNEL)
#@echo $(SLAW)
#@echo $(PP)
#@echo $(SC)
#@echo $(OBJS)

$(exec) : $(OBJS)
	$(FC) $(FCFLAGS) $^ -L$(NETCDF)/lib -lnetcdf -L$(NETCDFF)/lib -lnetcdff -o $@ 

%.o : %.f90
	$(FC) $(FCFLAG) -c -I$(NETCDFF)/include -I$(NETCDF)/include $< 

%.o : %.f
	$(FC) $(FCFLAG) -c $<

libcfgio.a : $(OBJS_IN)
	ar rc $@ $^; \
	ranlib $@; \

libspecial_functions.a :
	@echo 'libspecial_function build'
#	cd $(INC_SF); \
#	$(FC) -c f90split.f90; \
#	$(FC) f90split.o -o f90split; \
	mkdir temp; \
	cd temp; \
	../f90split ../special_functions.f90; \
	for FILE in `ls -1 *.f90`;\
 		do	$(FC) $(FCFLAGS) -c $$FILE; \
	done; \
	ar qc $@ *.o; \
	mv $@ ../../../; \
	cd ../; \
	rm -r ./temp/; \
	rm f90split.o; \
	cd ../../; \

clean:
	rm -f ./*.o ./*.a ./*.mod

cleano: 
	rm -f ./*.o ./*.mod