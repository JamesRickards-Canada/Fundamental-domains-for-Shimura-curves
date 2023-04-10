#Change this to where PARI/GP is installed.
PARI_LOC = /usr/local

#The files we want to include
SRCS = fdom.o fdom_extra.o
#fdom_extra.o Removing for now.

#Name of the output library
TARGET = fdom


#Nothing after here should need to be changed.

#More locations
PARI_LIB = $(PARI_LOC)/lib
PARI_INCLUDE = $(PARI_LOC)/include
PARI_CFG = $(PARI_LIB)/pari/pari.cfg

#Object names
OBJS = $(SRCS)
ALL = $(DYN)
VER = $(shell grep "pari_release=" "$(PARI_CFG)" | cut -d"'" -f2 | tr . - | cut -d"-" -f1,2)
DYN = lib$(TARGET)-$(VER).so

#Compiling options
CFLAGS     = -O3 -Wall -fno-strict-aliasing
CC         = cc
CPPFLAGS   = -I. -I$(PARI_INCLUDE)
LDFLAGS    = -O3 -Wall -fno-strict-aliasing    -Wl,--export-dynamic 
MODLDFLAGS = -shared  $(CFLAGS) $(DLCFLAGS) -Wl,-shared 
EXTRAMODLDFLAGS = -lc -lm -L$(PARI_LIB) -lpari
EXTRALIBS  =
DLCFLAGS   = -fPIC

#Recipes
all: $(DYN)

#Rule to build the library
$(DYN): $(OBJS)
	$(CC) -o $@ $(MODLDFLAGS) $(EXTRACFLAGS) $(OBJS) $(EXTRAMODLDFLAGS)

#Rule to make the object files
%.o: %.c
	$(CC) -c $(CFLAGS) $(EXTRACFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<

clean:
	-$(RM) *.o $(ALL)