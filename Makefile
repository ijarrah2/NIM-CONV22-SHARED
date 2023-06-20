# Fortran compiler
F90=gfortran
SILO_path=3rd_party/silo-4.11
# ============================================
F90SRC=src/main.f90 src/postprocess.f90
SRCOBJ=src/mesh.f90 src/parameters.f90 src/general.f90 src/restart.f90 src/curvilinear.f90 src/create_mesh.f90 src/create_geometry.f90 src/prepare_mesh.f90 src/geometry.f90 src/nim.f90 src/bc.f90 src/nse.f90
MODOBJ=$(SRCOBJ:.f90=.o)
PROGRAM=nim3d.exe
SILO_option=-I$(SILO_path)/include/ -L$(SILO_path)/lib/
OPTIONS= $(SILO_option) -lsilo -O3 -ffree-line-length-300 -fcheck=bounds -g -fopenmp -fdefault-real-8
# mkdir src/obj
# ============================================
$(PROGRAM):	$(MODOBJ)
	@echo $(MODOBJ)
	$(F90) $(OPTIONS) $(F90SRC) -lsilo -o $(PROGRAM) 
%.o:	%.f90
	$(F90) -c $(OPTIONS) $^ -o $@ 
clean:	
	rm -f $(MODOBJ) $(PROGRAM) *.mod plots.visit plots/step*
cleanall:
	rm -f $(MODOBJ) $(PROGRAM) *.mod plots.visit plots/step*
cleansilo:
	rm -r 3rd_party/silo-4.11
	rm 3rd_party/silo-4.11.tar.gz
