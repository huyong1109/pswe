FC= mpiifort
CPP = /usr/bin/cpp 
SRC= kinds_mod.f90 communicate.f90  pswe.f90
OBJS= $(SRC:.f90=.o) 
EXE=PSWE

#module_array.f90 cs.f90 dif.f90 euler.f90 haurwitz.f90 main.f90

all : $(SRC) $(EXE)

$(EXE): $(OBJS) 
	$(FC)  $(OBJS)  -o $@


#.SUFFIXES : .f90
.f90.o:
	$(FC) $<  -c $@
	#####> $*.f90 ; $(FC) -c $*.F90 ; mv $*.F90 $*.i
 

clean:
	rm -f *.o *.mod  $(EXE)
run:
	make 
	./$(EXE)
