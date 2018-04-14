#valgrind --leak-check=full --track-origins=yes ./prog
#CC			=icc
#CL			=icpc
#CFLAGS	= -O3 -g -fast -no-ipo -restrict -vec_report0 -ansi-alias -openmp -MMD -MP#-I ./
#CFLAGS	=-O3 -Wall -Werror-all -restrict -vec_report0 -parallel -openmp #-I ./
#CFLAGS	=-O0 -g -restrict -vec_report0 -parallel -openmp #-I ./

#CC		=g++
#CL		=g++
#CFLAGS 	=-O3 -msse -msse2 -msse3 -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=2 -fopenmp #-I ./
#CFLAGS 	= -O0 -g -D_ERROR_ -Wall -Weffc++ -fopenmp #-I ./
#CFLAGS		=-O3 -msse -msse2 -msse3 -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=0 -fopenmp -I ./
#LDFLAGS	= -lgsl -lgslcblas
#LDFLAGS	=

#CC       := g++
#CL       := g++
#CFLAGS 	 := -O3 -g -msse -msse2 -msse3 -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=0 -fopenmp -MMD -MP -finline-limit=20000 --param inline-unit-growth=1000 #-Winline
#CFLAGS 	 :=-g -fopenmp -MMD -MP
CC       := icc
CL       := icpc
CFLAGS 	 := -std=c++11 -O3 -g -no-ipo -restrict -fargument-noalias -ansi-alias -qopenmp -qopenmp-simd -MMD -MP #-opt-report=3#-Wall -Wextra -Winline -openmp-report0 -opt-report0 -vec_report3 -par_report3 -fast

VPATH    := efields:grids:models:observers
TARGET   := prog
LIBS     := 
EXT      := cpp
BUILDDIR := build
SRCDIRS  := efields grids models observers

#Petsc
ifeq ($(WITHPETSC), 1)
include ${PETSC_DIR}/lib/petsc/conf/variables
LIBS     += ${PETSC_LIB}
CFLAGS 	 += ${PETSC_CC_INCLUDES} -DPETSC
endif

override BUILDDIR := $(strip $(BUILDDIR))
override SRCDIRS  := $(addsuffix /MakefileAddSources,$(SRCDIRS))
SOURCES	 := main.cpp ConfigFile.cpp configuration.cpp system.cpp timeintegration.cpp
-include $(SRCDIRS)
#SOURCES := $(wildcard *.$(EXT))
OBJECTS  := $(addprefix $(BUILDDIR)/,$(SOURCES:.$(EXT)=.o))
DEPS     := $(addprefix $(BUILDDIR)/,$(SOURCES:.$(EXT)=.d))

$(TARGET): $(BUILDDIR) $(OBJECTS)
	$(CL) $(CFLAGS) -o $@ $(OBJECTS) $(LIBS)

$(BUILDDIR)/%.o : %.$(EXT)
	$(CC) $(CFLAGS) -o $@ -c $<

$(BUILDDIR):
	mkdir build

ifneq ($(MAKECMDGOALS), clean)
-include $(DEPS)
endif

clean:
	rm -rf $(BUILDDIR)
