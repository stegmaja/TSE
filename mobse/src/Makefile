IGNORE:
CMPLR = gfortran
FFLAGS = -O2 -fcheck=all
LFLAGS = ../input/const_mobse.h ../input/zdata.h 

.f.o:
	$(CMPLR) -c $(FFLAGS) $<
 
SRCE = \
comenv.f corerd.f deltat.f dgcore.f evolve.f gntage.f \
hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f \
star.f zcnsts.f zfuncs.f pisn.f eddington.f fallback.f

MAIN1 = mobse.f
OBJT1 = $(MAIN1:.f=.o) $(SRCE:.f=.o) 
mobse:	clean\
	$(OBJT1) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) -o mobse.x  

clean:
	rm -f *.o *.x
