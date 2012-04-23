#PROG = ~/bin/a.out  
 PROG = ./ciorca 

#--------------------------------------------------------------------------
 OSTYPE=LINUXI
#OSTYPE=LINUXL
#OSTYPE=LINUX    
# source /usr/qc/lf95/csh_setup nicht vergessen (Lahey compiler)
# source /usr/qc/intel/compiler70/ia32/bin/ifcvars.csh (Intel compiler)
#--------------------------------------------------------------------------

#OBJS1=freeze.o readl.o
# OBJS1=cefine.o susy.o
OBJS1=ciorca.o io.o

OBJS2 = 

OBJS = $(OBJS1) $(OBJS2)
#--------------------------------------------------------------------------


ifeq ($(OSTYPE),LINUXI)
  PREOPTS =
#  FC = gfortran  -static
  FC = ifort  -static -debug #-warn unused
  CC = gcc
#  LINKER = gfortran -static # -static-libcxa         
  LINKER = ifort -static # -static-libcxa         
#  LIBS    = -lmkl_ia32 -lguide -lpthread -L/usr/qc/intel/mkl/10.0/lib/32
# LIBS    = -lmkl_ia32 -lmkl_lapack -lmkl_solver -lguide -lpthread -L/usr/qc/intel/mkl701/lib/32
  PREFLAG = -E -P
  CCFLAGS = -O -DLINUX
  FFLAGS = -O #-w90 -O
endif                     

# diese ziele gibts:
.PHONY: all
.PHONY: clean
# dieses ist das erste auftretende,
# wird also beim aufruf von make erzeugt (default)
all: $(PROG)


#--------------------------------------------------------------------------
# example.f: printversion.h
# example.f haengt von printversion.h ab
# wenn sich also  printversion.h aendert, wird example.f
# (und damit auch example.o) neu gemacht.
# was auch geht:

#--------------------------------------------------------------------------
# implizite Regel zur Erzeugung von *.o aus *.F ausschalten
%.o: %.F


%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@
	
# aus *.F mache ein *.f
%.f: %.F
	@echo "making $@ from $<"
	$(CC) $(PREFLAG) $(PREOPTS) $< -o $@

# aus *.f mache ein *.o
%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

# aus *.c mache ein *.o
%.o: %.c
	@echo "making $@ from $<"
	$(CC) $(CCFLAGS) -c $< -o $@

# linken
$(PROG): $(OBJS) 
	$(LINKER) $(OBJS) $(LIBS) -o $(PROG)


#aufraeumen
clean:
	rm -f *.o $(PROG) 
#	rm -f $(patsubst %.F, %.f, $(wildcard *.F))

bak:
	cp io.f BAK/.
	cp ciorca.f BAK/.
	cp Makefile BAK/.
	cp *manual* BAK/.
sync:
	rsync -avzP ../CIADF v61u0036@lisa.sara.nl:
	rsync -auvzP ../CIADF kruse@lisa.sara.nl:
	rsync -auvzP ../CIADF hkruse@tc.few.vu.nl:
