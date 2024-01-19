#The standard location for the file pari.cfg is /usr/local/lib/pari/pari.cfg. If this is the case, call "make" to get going. Otherwise, call "make setup", which helps you find the correct pari.cfg file, saving it in paricfg_loc.txt (this step can be manually done as well, if you know where it is saved to). After this, call "make" to get going. Calling "make clean" clears up the .o files.

#The files we want to include
OBJS = fdom.o fdom_extra.o
#Name of the output library
TARGET = fdom

#Nothing after here should be modified, unless you know what you are doing.

#PARI_LIB is folder where libpari.so is found, PARI_INCLUDE is where the .h header files are found, and PARI_CFG is the location of pari.cfg.
PARI_LOC = paricfg_loc.txt
PARI_CFG = $(file < $(PARI_LOC))
ifeq ($(PARI_CFG), )
	PARI_CFG = /usr/local/lib/pari/pari.cfg
endif
PARI_LIB = $(shell grep "libdir=" "$(PARI_CFG)" -s | cut -d"'" -f2)
PARI_INCLUDE = $(shell grep "includedir=" "$(PARI_CFG)" -s | cut -d"'" -f2)/pari

#Naming the library file to include the version of pari/gp.
VER = $(shell grep "pari_release=" "$(PARI_CFG)" -s | cut -d"'" -f2 | tr . - | cut -d"-" -f1,2,3)
DYN = lib$(TARGET)-$(VER).so

#Compiling options
CC = cc
CFLAGS = -O3 -Wall -fno-strict-aliasing -fPIC
RM = rm -f

#Recipes
all: $(DYN)

#Build the shared library object
$(DYN): $(OBJS)
	$(CC) -o $@ -shared	$(CFLAGS) -Wl,-shared $(OBJS) -lc -lm -L$(PARI_LIB) -lpari

#Make the object files
%.o: %.c
	$(CC) -c $(CFLAGS) -I. -I$(PARI_INCLUDE) $<

#Clear all .o files
clean:
	$(RM) *.o $(ALL)

#Finds where pari/gp is installed and saves it to paricfg_loc.txt for the future.
setup:
	@printf "We need to find the location of the \e[33mconfiguration file\e[0m (\e[32mpari.cfg\e[0m), which contains the location of the \e[33mcompiled library\e[0m (\e[32mlibpari.so\e[0m) and the \e[33mheader files\e[0m. The standard location is:\n\t\e[32m/usr/local/lib/pari/pari.cfg\e[0m,\nbut this can depend on the system.\nWhere should we look for \e[0m \e[32mpari.cfg\e[0m? (default is \e[32m/usr\e[0m)"
	@read -r -p " " SEARCHLOC; \
	if [ "$$SEARCHLOC" = "" ] ; then \
		SEARCHLOC="/usr"; \
	fi; \
	printf "Searching for \e[32mpari.cfg\e[0m:\n"; \
	FOUNDLOCS=$$(find $$SEARCHLOC -name pari.cfg 2>/dev/null); \
	printf "Here are the places we found \e[32mpari.cfg\e[0m:\n\e[32m"; \
	for poss in $$FOUNDLOCS; do \
		printf "\t$$poss\n"; \
	done; \
	printf "\e[0mYou will now pick which one is correct.\n"; \
	CFGLOC=""; \
	for poss in $$FOUNDLOCS; do \
		printf "Is \e[35m$$poss\e[0m correct? "; \
		read -r -p "(y/n) " RESPONSE; \
		if [ "$$RESPONSE" = "y" ]; then \
			CFGLOC=$$poss; \
			break; \
		fi; \
	done; \
	if [ "$$CFGLOC" = "" ] ; then \
		printf "No file found or selected. Maybe PARI/GP is installed in a different directory than what you supplied?\n"; \
		exit; \
	fi; \
	printf "$${CFGLOC}" > ./paricfg_loc.txt; \
	echo "Setup complete!"
