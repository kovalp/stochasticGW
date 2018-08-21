## compiler options
FCMPI 	= mpif90
MPIFLG  = -DMPI -O3 
FFTFLG 	= -lfftw3


# compile
all :  clean sgw clnomod 

sgw: simmpi gwm atoms kbtop ppm kbmod
	$(FCMPI) $(MPIFLG)  -o sgw.x kb/*f90 libgw/*f libgw/*f90 *f90  $(FFTFLG) 

simmpi:
	$(FCMPI) $(MPIFLG) -c ./libgw/0_library_mpi_module.f90 ./libgw/0_simple_mpi.f90

gwm:
	$(FCMPI) $(MPIFLG) -c 1_gwm.f90

atoms:
	$(FCMPI) $(MPIFLG) -c 2_atoms.f90

kbtop:
	$(FCMPI) $(MPIFLG) -c 2_kb_top_module.f90

ppm:
	$(FCMPI) $(MPIFLG) -c 3_ppm.f90

kbmod:
	$(FCMPI) $(MPIFLG) -c kb/1_kb_mod.f90

clnomod:
	rm -rf *o *mod kb/*o libgw/*o 

clean:
	rm -rf *o *mod kb/*o libgw/*o sgw.x 
