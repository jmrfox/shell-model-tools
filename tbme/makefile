CFLAGS   = -O   # -CB
LFLAGS   =  # -CB 
DPLIBS   = 
COMPILER = ftn

TBME_OBJECTS = \
tbme_modules.o \
makeTBME.o \
subTBME.o \
interactions.o \
sub_contact_integral.o \
cwjhofcns.o \
talmi.o \
subHcm.o \
libra.o

%.o : %.f
	$(COMPILER) $(CFLAGS) -c $(<F)
%.o : %.f90
	$(COMPILER) $(CFLAGS) -c $(<F)



tbme.x : $(TBME_OBJECTS)
	$(COMPILER) -o tbme.x $(TBME_OBJECTS) $(LFLAGS) $(DPLIBS)
clean :
	rm tbme.x $(TBME_OBJECTS)
