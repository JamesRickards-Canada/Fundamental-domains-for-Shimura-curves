#Start with some constant definitions
PARI_LIB = /usr/local/lib
PARI_INCLUDE = /usr/local/include

SRCS = fdom.o fdom_extra.o
OBJS = $(SRCS)

TARGET = fdom
SHELL  = /bin/sh

CFLAGS     = -O3 -Wall -fno-strict-aliasing
EXTRACFLAGS=

CC         = cc
CPPFLAGS   = -I. -I$(PARI_INCLUDE)
LD         = cc
LDFLAGS    = -O3 -Wall -fno-strict-aliasing    -Wl,--export-dynamic 
MODLD      = cc
MODLDFLAGS = -shared  $(CFLAGS) $(DLCFLAGS) -Wl,-shared 
EXTRAMODLDFLAGS = -lc -lm -L$(PARI_LIB) -lpari
EXTRALIBS  =

DLCFLAGS   = -fPIC

DYN = lib$(TARGET).so
ALL = $(DYN)

all: $(DYN)

#Rule to build the library
$(DYN): $(OBJS)
	$(MODLD) -o $@ $(MODLDFLAGS) $(EXTRACFLAGS) $(OBJS) $(EXTRAMODLDFLAGS)

#Rule to make the object files
%.o: %.c
	$(CC) -c $(CFLAGS) $(EXTRACFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<

clean:
	-$(RM) *.o $(ALL)
