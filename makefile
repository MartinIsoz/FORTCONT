#=====COMPILER==========================================================
f90comp=gfortran
#f77comp=gfortran
#=====FLAGS=============================================================
# flags for debugging or maximal performace
f90compFLAGS= -g -Wall -fbounds-check
#f90compFLAGS= -O2
# flags for all - look for .mod files
f90compFLAGS+= -I/home/martin/Documents/Phd/01_BADS/GNUFOR2 
#library flags
libFLAGS 	= -L/usr/local/lib 
libs		= -llapack -lblas
#=====PROGRAM=NAMES=====================================================
main		= testfortcont
#=====CREATE=RESULT=PROGRAM=============================================
$(main).out: $(main).o vecfields.o fortcont.o gause.o
	$(f90comp)  $(f90compFLAGS) $(main).o vecfields.o fortcont.o gause.o $(libFLAGS) $(libs) -o testfortcont.out
#=====CREATE=OBJECTS====================================================
fortcont.o: fortcont.f90
	$(f90comp) -c $(f90compFLAGS) fortcont.f90 $<
vecfields.o:vecfields.f90
	$(f90comp) -c $(f90compFLAGS) vecfields.f90 $<
$(main).o:  $(main).f90 fortcont.o vecfields.o
	$(f90comp) -c $(f90compFLAGS) $(main).f90 $<
#=======================================================================
clean:
	rm *.o *.mod
#=======================================================================
run: 
	./$(main).out
#=======================================================================
