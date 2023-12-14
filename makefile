
FC	= gfortran
CC	= gcc
LD	= $(FC)
F77	= $(FC)
SWP	= swplist
RM	= /bin/rm -f
MAKE    = smake -J 8
MP	=
ABI	= 
ISA	= 
PROC	= 
ARCH	= $(MP) $(ABI) $(ISA)
OLEVEL  =  -O3
#OLEVEL  =  -O3 -scalar_rep  -vec -prefetch -tpp6
#OLEVEL  =  -O3 -scalar_rep  -vec -prefetch -tpp6  -axK -pad
#OLEVEL	= -g -DEBUG:conform_check=ON:div_check=3:subscript_check=ON:trap_uninitialized=ON:varargs_interface_check=ON:verbose_runtime=ON
DEBUGFLAGS = "-g -DEBUG:conform_check=ON:div_check=3:subscript_check=ON:trap_uninitialized=ON:varargs_interface_check=ON:verbose_runtime=ON"

FOPTS   = -r8  

COPTS	= 
FFLAGS	= $(ARCH) $(OLEVEL) $(FOPTS)
CFLAGS	= -g $(ARCH) $(COPTS)
#LIBS   = -lPEPCF90
LIBS    =
LDFLAGS	= $(ARCH) $(OLEVEL)
PROF	=

FOBJS   = dims.o mesh_var.o flo_var.o iodefs.o\
	  time_var.o flo_param.o solv_var.o\
	  alloc_mesh.o alloc_flo.o\
	  input.o mesh.o init.o tstep.o\
	  update.o bconds.o\
	   derivs_bgk.o \
	  bgkflux.o \
	  main.o output.o


COBJS	= 
OBJS	= $(FOBJS) $(COBJS)
EXEC	= mlayer

$(EXEC):	$(OBJS)
	$(LD) -o $@ $(PROF) $(LDFLAGS) $(QV_LOPT) $(OBJS) $(LIBS)


debug:
	$(MAKE) -e OLEVEL=$(DEBUGFLAGS) $(EXEC)

main.o: main.f90
	$(FC)  $(CFLAGS) -c main.f90

dims.o: dims.f90
	$(FC)  $(CFLAGS) -c dims.f90

time_var.o: time_var.f90
	$(FC)  $(CFLAGS) -c time_var.f90

iodefs.o: iodefs.f90
	$(FC)  $(CFLAGS) -c iodefs.f90

mesh_var.o: mesh_var.f90
	$(FC)  $(CFLAGS) -c mesh_var.f90

flo_param.o: flo_param.f90
	$(FC)  $(CFLAGS) -c flo_param.f90

input.o: input.f90
	$(FC)  $(CFLAGS) -c input.f90

flo_var.o: flo_var.f90
	$(FC)  $(CFLAGS) -c flo_var.f90

mesh.o: mesh.f90
	$(FC)  $(CFLAGS) -c mesh.f90

alloc_mesh.o: alloc_mesh.f90
	$(FC)  $(CFLAGS) -c alloc_mesh.f90

init.o: init.f90
	$(FC)  $(CFLAGS) -c init.f90

alloc_flo.o: alloc_flo.f90
	$(FC)  $(CFLAGS) -c alloc_flo.f90

solv_var.o: solv_var.f90
	$(FC)  $(CFLAGS) -c solv_var.f90

bconds.o: bconds.f90
	$(FC)  $(CFLAGS) -c bconds.f90

output.o: output.f90
	$(FC)  $(CFLAGS) -c output.f90

tstep.o: tstep.f90
	$(FC)  $(CFLAGS) -c tstep.f90

update.o: update.f90
	$(FC)  $(CFLAGS) -c update.f90

derivs_bgk.o: derivs_bgk.f90
	$(FC)  $(CFLAGS) -c derivs_bgk.f90

bgkflux.o: bgkflux.f90
	$(FC)  $(CFLAGS) -c bgkflux.f90



clean:
	$(RM) $(EXEC) $(OBJS)
	rm -f *.log *~ *.dat *.mod *.rst *.tec fort.*

#.SUFFIXES:
.SUFFIXES: .o .F .c .f *.f90 .swp

.F.o:
	$(FC)  -c $(FFLAGS) $(QV_OPT) $(DEFINES) $<

.f90.o:
	$(FC)  -c $(FFLAGS) $(QV_OPT) $(DEFINES) $<

.f.o:
	$(FC)  -c $(FFLAGS) $(QV_OPT) $(DEFINES) $<

.c.o:
	$(CC)  -c $(CFLAGS) $(QV_OPT) $(DEFINES) $<

.F.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.f.swp:
	$(SWP) -c $(FFLAGS) $(QV_OPT) $(DEFINES) -WK,-cmp=$*.m $<

.c.swp:
	$(SWP) -c $(CFLAGS) $(QV_OPT) $(DEFINES) $<
