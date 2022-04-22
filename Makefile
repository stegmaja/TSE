CMPLR = gfortran 
FFLAGS = -O2 $(profile) $(GCDA) $(coverage)
FFLAGS = -O2 
LFLAGS = const_bse.h zdata.h 

.f.o:
	$(CMPLR) -c $(FFLAGS) $<

SRCE1 = common.f sub.f supernova.f evolve.f deltat.f evolv1.f hrdiag.f kick.f mlwind.f mrenv.f ran3.f star.f zcnsts.f zfuncs.f comenv.f corerd.f gntage.f instar.f rl.f dgcore.f mix.f evolv2.f bse.f sse.f sse-rl.f LKspin-s.f
SRCE2 = common.f sub.f supernova.f evolve.f deltat.f evolv1.f hrdiag.f kick.f mlwind.f mrenv.f ran3.f star.f zcnsts.f zfuncs.f comenv.f corerd.f gntage.f instar.f rl.f dgcore.f mix.f evolv2.f bse.f sse.f sse-rl.f initialise.f

OBJT1 = $(SRCE1:.f=.o)
OBJT2 = $(SRCE2:.f=.o)

tripleStars: $(OBJT1) $(LFLAGS)
	$(CMPLR)   $(FFLAGS) $(OBJT1) -o LKspin-s

initialise: $(OBJT2) $(LFLAGS)
	$(CMPLR)   $(FFLAGS) $(OBJT2) -o initialise

clean:
	rm -rf LKspin-s initialise *.o *.mod ./0.*/ *.gcno
