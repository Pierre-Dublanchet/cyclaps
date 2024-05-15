#--cyclaps makefile---#
FC=mpif90
FCFLAGS=-O5
#--Indicate path to netcdf librairies below--#
NETCDF=/usr
NETCDFF=/usr
#--------------------------------------------#
INC_SF=./src/special_functions/
GEOM=$(word 2, $(subst _, ,$(exec)))
KERNEL=$(word 3, $(subst _, ,$(exec)))
SLAW=$(word 4, $(subst _, ,$(exec)))
PP=$(word 5, $(subst _, ,$(exec)))
SC=$(word 6, $(subst _, ,$(exec)))

SRC_DIRS = $(sort $(dir ./src/ ./src/cfgio/ ./src/fft/ ./src/special_functions/))
vpath %.f90 $(SRC_DIRS)
vpath %.f $(SRC_DIRS)

ifeq "$(KERNEL)" "freesurface"
	OBJS := mod_cst_fault_rns.o string_conv_mod.o cfgio_mod.o
else
	OBJS := mod_cst_fault_rns.o string_conv_mod.o cfgio_mod.o zfft1d.o factor.o pfactor.o fft235.o kernel.o
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

ifeq "$(KERNEL)" "cr"
	OBJS :=  $(OBJS) libspecial_functions.a
endif

do: $(exec)

$(exec) : $(OBJS)
	#$(FC) $(FCFLAGS) $^ -L$(NETCDF)/lib -lnetcdf -L$(NETCDFF)/lib -lnetcdff -o $@ 
	$(FC) $(FCFLAGS) $^ -lnetcdf -lnetcdff -o $@ 

%.o : %.f90
	$(FC) $(FCFLAG) -c -I$(NETCDFF)/include -I$(NETCDF)/include $< 

%.o : %.f
	$(FC) $(FCFLAG) -c $<

libcfgio.a : $(OBJS_IN)
	ar rc $@ $^; \
	ranlib $@; \

libspecial_functions.a :
	@echo 'libspecial_function build'
	cd $(INC_SF); \
	$(FC) -c f90split.f90; \
	$(FC) f90split.o -o f90split; \
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
