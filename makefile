#***********************************************************************
# Makefile for PHITS, last revised 2014/08/29
#***********************************************************************

### Machine Dependent variables, please set your environment 
ENVFLAGS = MacGfort
# Linux In1tel Fortran: LinIfort
# Windows gfortran (4.71 or later, OpenMP unstable): WinGfort
# Windows Cygwin gfortran (4.7 or later, OpenMP unstable): WinCygGfort
# Mac Intel fortran: MacIfort
# Mac gfortran (4.71 or later): MacGfort
# BX900 (Fujitsu compiler): BX900
# FX1 (you have to use "gmake" instead of "make"): FX1
### If you succeeded in compiling PHITS in new environment, 
### please send your makefile to phits-office@jaea.go.jp

# If you want to use MPI, delete # in the next line
# MPIFLAGS = true
# If you want to use OpenMP, delete # in the next line
# OMPFLAGS = true

# If you use option '-j', you will get many errors because the order of compiling
# files should not be changed in the make of PHITS. ('-j' option changes the order)
# In that case, you have to type 'make -j' again, or do not use '-j' option
# then you will succeed in making the PHITS executable.

ifeq ($(ENVFLAGS),)
$(error Error! ENVFLAGS is not defined)
endif

### Linux Intel Fortran
ifeq ($(ENVFLAGS),LinIfort)
 TARGET = ../phits_lin.exe
 FC      = mpif90
 ifeq ($(OMPFLAGS),true) 
  OMPFLAGS = -openmp
 endif
 ifeq ($(MPIFLAGS),true)
  FC = mpif90
  SRCS8   = mpi-lin.f unix90.f mdp-uni90.f
  LIBS = -lmpi $(OMPFLAGS)  ### -lmpi sometimes -lmpich, or not necessary
 else
  FC  = ifort
  CC  = icc
  CXX = icpc
  SRCS8   = mpi-non.f unix90.f mdp-uni90.f
  LIBS = $(OMPFLAGS) -lstdc++
 endif
 INCLUDES =
 FCFLAGS = -O3 -fpconstant -fno-second-underscore
 CXXFLAGS= -03 `root-config --cflags`
 CXXLIBS  = `root-config --libs`
 CCFLAGS  = -O3
 CCLIBS   =
 LD      = $(FC) 
 LDFLAGS =
endif

### BX900
ifeq ($(ENVFLAGS),BX900)
limit vmemoryuse 2097152 # extend the limitation of memory usage
#ulimit -v 2097152        # for bash user, replace above command by this line
TARGET  = ../phits_bx900.exe
ifeq ($(OMPFLAGS),true) 
 OMPFLAGS = -KOMP
endif
ifeq ($(MPIFLAGS),true)
 FC      = mpifrt
 SRCS8   = mpi-lin.f unix90.f mdp-uni.f
else
 FC      = frt
 SRCS8   = mpi-non.f unix90.f mdp-uni.f
endif
 FCFLAGS = -fw -X9 -Am -Kfast,ocl,preex,noomitfp $(OMPFLAGS) \
           -xzdistfunc,Gau -xtime_b -xshellspar -xfnf \
           -xsmmmod.xi,smmmod.eps
LIBS    =
LD      = $(FC)
LDFLAGS = $(OMPFLAGS)
endif

### FX1, OpenMP uncheck, you have to use "gmake" instead of "make"
ifeq ($(ENVFLAGS),FX1)
TARGET  = ../phits_fx1.exe
ifeq ($(OMPFLAGS),true) 
 OMPFLAGS = -KOMP
endif
ifeq ($(MPIFLAGS),true)
 FC      = mpifrt
 SRCS8   = mpi-lin.f unix90.f mdp-uni.f
else
 FC      = frt
 SRCS8   = mpi-non.f unix90.f mdp-uni.f
endif
OPTFLAGS = -Kfast,ocl,preex
FCFLAGS  = -fw -X9 -Am $(OMPFLAGS) \
           -xzdistfunc,Gau -xtime_b -xshellspar -xfnf \
           -xsmmmod.xi,smmmod.eps
LIBS    =
LD      = $(FC)
LDFLAGS = $(OMPFLAGS)
endif

### Mac Intel fortran
ifeq ($(ENVFLAGS),MacIfort)
TARGET   = ../phits_mac.exe
FC       = ifort
ifeq ($(OMPFLAGS),true)
  OMPFLAGS  = -openmp
  LD = $(FC) -static-intel
else
  OMPFLAGS  = 
  LD = $(FC)
endif
ifeq ($(MPIFLAGS),true)
  FCFLAGS  = $(OMPFLAGS) -fpconstant -O3
  LDFLAGS  = $(OMPFLAGS) -Wl,-commons,use_dylibs -I/usr/local/lib
  INCLUDES = -I/usr/local/include -Wl,-commons,use_dylibs -I/usr/local/lib
  LIBS     = -L/usr/local/lib -lmpi_f90 -lmpi_f77 -lmpi -lm
  SRCS8    = mpi-lin.f unix90.f mdp-uni90.f
else
  FCFLAGS  = $(OMPFLAGS) -fpconstant -O3 `root-config --cflags`
  LDFLAGS  = $(OMPFLAGS) `root-config --ldflags --libs`
  INCLUDES = #-I/Users/youkanau/Documents/root/include
  LIBS     = #-L/Users/youkanau/Documents/root/lib
  SRCS8    = mpi-non.f unix90.f mdp-uni90.f
endif
endif

### Mac gfortran (4.71 or later)
ifeq ($(ENVFLAGS),MacGfort)
TARGET   = ../phits270.exe
FC       = gfortran-mp-4.8
CC       = gcc
CXX      = clang++
LD       = $(FC)
ifeq ($(OMPFLAGS),true)
  OMPFLAGS = -fopenmp
else
  OMPFLAGS =
endif
ifeq ($(MPIFLAGS),true)
  FCFLAGS  = $(OMPFLAGS) -O0 $(shell mpif90-openmpi-mp -showme:compile) -fdefault-double-8 -fdefault-real-8 
  LDFLAGS  = $(OMPFLAGS) $(shell mpif90-openmpi-mp -showme:link)
  INCLUDES = -I$(shell mpif90-openmpi-mp -showme:incdirs)
  LIBS     = -lmpi_usempi -lmpi_mpifh -lmpi
  SRCS8    = mpi-lin.f unix90.f mdp-uni90.f
else
  #VPATH    = usrsrc
  FCFLAGS  = $(OMPFLAGS) -O0 -fdefault-double-8 -fdefault-real-8 -c
  LDFLAGS  = $(OMPFLAGS) -lstdc++ #-I include
  CCFLAGS  = -O3   
  CXXFLAGS = -g `root-config --cflags`
  INCLUDES = 
  LIBS =
  CXXLIBS     = -L$(ROOTSYS)/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl
  SRCS8    = mpi-non.f unix90.f mdp-uni90.f
endif
endif

### Windows gfortran (4.71 or later, MPI unchecked, OpenMP unstable)
### You have to install "gfortran".
### Download site for GCC including gfortran for Windows
### (http://tdm-gcc.tdragon.net/download) Bundle installer is recommended.
### This package also has "mingw32-make" you can use make command for Windows.
ifeq ($(ENVFLAGS),WinGfort)
TARGET   = ../phits_win.exe
FC       = gfortran
LD       = $(FC)
ifeq ($(OMPFLAGS),true)
  OMPFLAGS = -fopenmp
else
  OMPFLAGS =
endif
ifeq ($(MPIFLAGS),true)
  FCFLAGS  = $(OMPFLAGS) -O1 -fdefault-double-8 -fdefault-real-8
  # LDFLAGS  = $(OMPFLAGS) -Wl,-commons,use_dylibs -I/usr/local/lib
  # INCLUDES = -I/usr/local/include -Wl,-commons,use_dylibs -I/usr/local/lib
  # LIBS     = -L/usr/local/lib -lmpi_f90 -lmpi_f77 -lmpi -lm
  SRCS8    = mpi-lin.f unix90.f mdp-win.f
else
  FCFLAGS  = $(OMPFLAGS) -O1 -fdefault-double-8 -fdefault-real-8
  LDFLAGS  = $(OMPFLAGS)
  INCLUDES =
  LIBS     =
  SRCS8    = mpi-non.f unix90.f mdp-win.f
endif
endif

### Cygwin gfortran (4.7 or later, OpenMP unstable, MPI unchecked)
ifeq ($(ENVFLAGS),WinCygGfort)
TARGET   = ../phits_cyg.exe
FC       = gfortran
LD       = $(FC)
ifeq ($(OMPFLAGS),true)
  OMPFLAGS = -fopenmp
else
  OMPFLAGS =
endif
ifeq ($(MPIFLAGS),true)
  FCFLAGS  = $(OMPFLAGS) -O0 -fdefault-double-8 -fdefault-real-8
  # LDFLAGS  = $(OMPFLAGS) -Wl,-commons,use_dylibs -I/usr/local/lib
  # INCLUDES = -I/usr/local/include -Wl,-commons,use_dylibs -I/usr/local/lib
  # LIBS     = -L/usr/local/lib -lmpi_f90 -lmpi_f77 -lmpi -lm
  SRCS8    = mpi-lin.f unix90.f mdp-win.f
else
  FCFLAGS  = $(OMPFLAGS) -O0 -fdefault-double-8 -fdefault-real-8
  LDFLAGS  = $(OMPFLAGS) -Wl,--stack,8388608
  INCLUDES =
  LIBS     =
  SRCS8    = mpi-non.f unix90.f mdp-win.f
endif
endif


#=======================================================================
# modules
#=======================================================================
SRCS0 = talmod.f \
				mmbankmod.f membankmod.f ggbankmod.f \
				ggmbankmod.f eventtalmod.f restalmod.f \
				smm.f gammod2.f gammod1.f gammod.f ggm05mod.f \
				gammod-a.f gammod-b.f gammod-c.f gammod-d.f gammod-e.f \
				egs5mod.f 

#=======================================================================
# for machine dependent or user defined source and analysis
#=======================================================================
SRCS1 = \
				usrsors.f  usrmgf1.f  usrmgf3.f  anal-002.f \
				usrdfn1.f  usrdfn2.f


#=======================================================================
# for param.inc
#=======================================================================
SRCS2 = \
				analyz.f   celimp.f   dataup.f   getflt.f   magtrs.f   \
				nreac.f    ovly12.f   ovly13.f   partrs.f   range.f    \
				read00.f   read01.f   read02.f   sors.f     talls00.f  \
				talls01.f  talls02.f  talls03.f  talls04.f  talls05.f  \
				talls06.f  talls07.f  tallsm1.f  tallsm2.f  tallsm3.f  \
				update.f   wrnt12.f   wrnt13.f   read03.f   marscg.f   \
				ggs00.f    ggs01.f    ggs02.f    ggs03.f    wrnt10.f   \
				geocntl.f  ggm01.f    ggm02.f    ggm03.f    ggm04.f    \
				ggm05.f    ggm06.f    ggm07.f    ggm08.f    a-angel.f  \
				ovly14.f   ovly15.f   usrelst1.f usrelst2.f usrmgt1.f  \
				usrmgt2.f  talls08.f  statistic.f tallsm4.f itrminmax.f

#=======================================================================
# for new package
#=======================================================================
SRCS3 = \
				main.f     dklos.f   ncasc.f     nelst.f    nevap.f    \
				sdml.f     gem.f     gemset.f    utl01.f    \
				utl02.f    jbook.f   masdis.f    atima01.f  atima02.f  \
				atima03.f  fismul.f  kurotama0.f usrtally.f incelf.f   \
				incelfin.f inclp.f   inclin.f    dwbain.f   dwbaD.f    \
				dwbaE.f    photnucl.f

#=======================================================================
# for old package
#=======================================================================
SRCS4 = \
				bert.f     bertin.f   bert-bl0.f bert-bl1.f bert-bl2.f \
				utlnmtc.f  gamlib.f   erupin.f   erup.f     fissn.f    \
				isobert.f  isodat.f   randmc.f   energy.f   ndata01.f  \
				mars00.f   mars01.f   mars02.f   mars03.f   mars04.f

#=======================================================================
# for JAM
#=======================================================================
SRCS5 = \
				jamin.f    jam.f      jamdat.f   jamcoll.f  jamdec.f   \
				jamcross.f jampdf.f   jamsoft.f  jamhij.f   jamhard.f  \
				jambuu.f   jamana.f   pyjet.f    pythia.f   pysigh.f

#=======================================================================
# for JQMD
#=======================================================================
SRCS6 = \
				qmd00.f    qmdcoll.f  qmddflt.f  qmdgrnd.f  qmdinit.f   \
				qmdmfld.f                                               \
				qmd00_qmd-nr.f     qmdcoll_qmd-nr.f    qmdgrnd_qmd-nr.f \
				qmdmfld_qmd-nr.f

#=======================================================================
# for ANGEL except for utl03.f, a-func.f, a-utl00.f
#=======================================================================
SRCS7 = \
				utl03.f    a-func.f   a-utl00.f     \
				a-main0.f  a-main1.f  a-hsect.f  a-line.f   a-wtext.f

#=======================================================================
# for restart function
#=======================================================================
SRCS9 = \
				restart.f resutl.f \
				resttrack.f restcross.f restheat.f restdeposit.f restyield.f \
				restlet.f   restsed.f   restdpa.f  restproduct.f resttime.f  \
				restdeposit2.f reststar.f 

#======================================================================= 
# for EGS5
#======================================================================= 
SRCS10 = \
				egs5.o  egs5init.o      egs5pfpl.o      egs5pcoll.o     \
				egs5efpl.o      egs5ecoll.o     egs5ede.o       egs5annih.o     \
				egs2phits.o     egs5esteps0.o	egs5edxde.o	egs5init2.o	\
				egs5elast.o     egs5ededx.o

#=======================================================================
# for C file
#=======================================================================
SRCS11 = \
 				CoilMapReader.c

#=======================================================================
# for C++ file
#=======================================================================
SRCS12 = \
				test.cpp

OBJS0 = $(SRCS0:.f=.o)
OBJS1 = $(SRCS1:.f=.o)
OBJS2 = $(SRCS2:.f=.o)
OBJS3 = $(SRCS3:.f=.o)
OBJS4 = $(SRCS4:.f=.o)
OBJS5 = $(SRCS5:.f=.o)
OBJS6 = $(SRCS6:.f=.o)
OBJS7 = $(SRCS7:.f=.o)
OBJS8 = $(SRCS8:.f=.o)
OBJS9 = $(SRCS9:.f=.o)
OBJS10 = $(SRCS10:.f=.o)
OBJS11 = $(SRCS11:.c=.o)
OBJS12 = $(SRCS12:.cpp=.o)

OBJS = $(OBJS0) $(OBJS1) $(OBJS2) $(OBJS3) \
			 $(OBJS4) $(OBJS5) $(OBJS6) \
			 $(OBJS7) $(OBJS8) $(OBJS9) \
			 $(OBJS10) $(OBJS11) $(OBJS12)

$(TARGET) : $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $^ $(CXXLIBS) $(LIBS) 

.f.o:
	$(FC) -c $(FCFLAGS) $(OPTFLAGS) $*.f 
.c.o:
	$(CC) -c $(CCFLAGS) $*.c
.cpp.o:
	$(CXX) -c $(CXXFLAGS) $*.cpp
gammod.o: gammod.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
gammod1.o: gammod1.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
gammod-a.o: gammod-a.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
gammod-b.o: gammod-b.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
gammod-c.o: gammod-c.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
gammod-d.o: gammod-d.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
gammod-e.o: gammod-e.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)

ifeq ($(ENVFLAGS),WinGfort)
  read00.o: read00.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
  read03.o: read03.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
endif
ifeq ($(ENVFLAGS),WinCygGfort)
  read00.o: read00.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
  read03.o: read03.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
endif
ifeq ($(ENVFLAGS),MacGfort)
  read00.o: read00.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
  read03.o: read03.f
	$(FC) -c $(FCFLAGS) -O0 $(INCLUDES) $*.f $(LIBS)
endif

$(OBJS2) : param.inc
$(OBJS7) : angel01.inc
$(OBJS8) : param.inc

clean:
ifeq ($(ENVFLAGS),WinGfort)
	@del -f $(OBJS) *.mod
else
	@rm -f $(OBJS) *.mod
endif
