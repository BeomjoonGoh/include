# C++ compiler
CC = /opt/rh/devtoolset-2/root/usr/bin/gcc
#CC = /opt/intel/Compiler/13.1/composer_xe_2013.5.192/bin/intel64/icc

# complier flags
CF_CPP		= -std=c++11 
CF_WAR		= -Wall -Wno-non-template-friend -Wno-sign-compare 
CF_OMP		= -fopenmp 
CF_DEF		= -D NDEBUG
CF_OPT		= -O3
#CF_OPT		= -O1
#CF_DEBUG 	= -p -g

CFLAGS 		= $(CF_CPP) $(CF_WAR) $(CF_OMP) $(CF_OPT)# $(CF_DEF) 
#CFLAGS 	= $(CF_CPP) $(CF_WAR) $(CF_OPT) $(CF_DEBUG) 

# includes
INC_MYCPP 	= -I/home/bjgoh1990/lib/cpp/
INC_SPECIAL 	= -I/home/bjgoh1990/lib/cpp/specialfunctions/
INC_CBLAS 	= -I/home/bjgoh1990/lib/lapack-3.8.0/CBLAS/include/
INC_LAPACKE 	= -I/home/bjgoh1990/lib/lapack-3.8.0/LAPACKE/include/ 

INCLUDES = $(INC_MYCPP) $(INC_SPECIAL) $(INC_CBLAS) $(INC_LAPACKE) 

# library paths
LF_MYCPP	= -L/home/bjgoh1990/lib/cpp
LF_LAPACK	= -L/home/bjgoh1990/lib/lapack-3.8.0

LFLAGS 		= $(LF_LAPACK)

# libraries to link (the order matters!) libfoo.a or libfoo.so => -lfoo
LIB_CPP		= -lstdc++ -lm
LIB_OMP		= -lgomp
LIB_CBLAS	= -lcblas 
LIB_LAPACK 	= -llapacke -llapack -lrefblas -ltmglib -lgfortran
#LIB_MKL	= -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
#LIB_PT		= -lpthread

LIBS 		= $(LIB_CPP) $(LIB_OMP) $(LIB_CBLAS) $(LIB_LAPACK) 
#LIBS 		= $(LIB_CPP) $(LIB_OMP) $(LIB_MKL) $(LIB_PT)

# source and objects
CPPFILE = main
SRCS = $(CPPFILE).cpp
OBJS = $(SRCS:.cpp = .o)

# the executable file
MAIN = $(CPPFILE)
#MAIN = out

## Makefile ##
all: 	$(MAIN)

%.o: %.cpp  $(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(MAIN): $(OBJS)
	 $(CC) $(CFLAGS) $(INCLUDES) \
	 -o $@ $^ $(LFLAGS) $(LIBS)


RM = rm
clean:
		$(RM) -f $(MAIN)


## This is a suffix replacement rule for building .o's from .c's
## it uses automatic variables $<: the name of the prerequisite of
## the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
## (see the gnu make manual section about automatic variables)
#.c.o:
#		$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@
#
#clean:
#		$(RM) *.o *~ $(MAIN)
#
#depend: $(SRCS)
#		makedepend $(INCLUDES) $^
#
## DO NOT DELETE THIS LINE -- make depend needs it
