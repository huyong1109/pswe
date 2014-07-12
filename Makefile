FC= mpiifort 
CPP = /usr/bin/cpp 
SRC= kinds_mod.f90 communicate.f90  \
     exit_mod.f90 broadcast.f90\
     module_para.f90 distribution.f90 \
     module_array.f90 global_reductions.f90 \
     boundary.f90 module_io.f90\
     cs.f90 haurwitz.f90 \
     dif.f90 euler.f90 \
     pswe.f90 
OBJS= $(SRC:.f90=.o) 
EXE=PSWE


all : $(SRC) $(EXE)

$(EXE): $(OBJS) 
	$(FC)  $(OBJS)  -o $@


#.SUFFIXES : .f90
%.o: %.f90
	$(FC) -c -g -traceback   $<
 

clean:
	rm -f *.o *.mod  $(EXE)
run:
	mpirun -np 20 ./$(EXE)
