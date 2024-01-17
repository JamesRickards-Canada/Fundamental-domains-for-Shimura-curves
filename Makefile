#If the PARI/GP version you want is not installed in /usr/local, call "make setup" to create pari_loc.txt, which helps you find the correct installation locations. After this, call make to get going. make clean clears up the files.

#The files we want to include
OBJS = fdom.o fdom_extra.o
#Name of the output library
TARGET = fdom
#Nothing after here should be modified, unless you know what you are doing.

#PARI_LIB is where libpari.so is found, PARI_INCLUDE is where the .h header files are found, and PARI_CFG is where pari.cfg is found.
PARI_ALT_LOC = pari_loc.txt
ALT_LOCK_PATHS = $(file < $(PARI_ALT_LOC))
ifeq ($(ALT_LOCK_PATHS), )
	PARI_LIB = /usr/local/lib
	PARI_INCLUDE = /usr/local/include/pari
else
	PARI_LIB = $(word 1, $(ALT_LOCK_PATHS))
	PARI_INCLUDE = $(word 2, $(ALT_LOCK_PATHS))
endif
PARI_CFG = $(PARI_LIB)/pari/pari.cfg

#Naming the library file to include the version of pari/gp.
VER = $(shell grep "pari_release=" "$(PARI_CFG)" -s | cut -d"'" -f2 | tr . - | cut -d"-" -f1,2)
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

#Finds where pari/gp is installed and saves it to pari_loc.txt for the future.
setup:
	@printf "We need to find the location of the \e[33mcompiled library\e[0m (\e[32mlibpari.so\e[0m), the \e[33mconfiguration file\e[0m (\e[32mpari.cfg\e[0m), and the \e[33mheader files\e[0m. The standard locations are in:\n\t\e[32m/usr/local/lib/\n\t/usr/local/include\e[0m,\nbut this can depend on the system.\nWhere should we look for the files? (default is \e[32m/usr\e[0m)"
	@read -r -p " " SEARCHLOC; \
	if [ "$$SEARCHLOC" = "" ] ; then \
		SEARCHLOC="/usr"; \
	fi; \
	printf "Searching for \e[32mlibpari.so\e[0m:\n"; \
	FOUNDLOCS=$$(find $$SEARCHLOC -name libpari.so 2>/dev/null); \
	printf "Here are the places we found \e[32mlibpari.so\e[0m:\n\e[32m"; \
	for poss in $$FOUNDLOCS; do \
		printf "\t$$poss\n"; \
	done; \
	printf "\e[0mYou will now pick which one is correct.\n"; \
	LIBLOC=""; \
	for poss in $$FOUNDLOCS; do \
		printf "Is \e[35m$$poss\e[0m correct? "; \
		read -r -p "(y/n) " RESPONSE; \
		if [ "$$RESPONSE" = "y" ]; then \
			LIBLOC=$$poss; \
			break; \
		fi; \
	done; \
	if [ "$$LIBLOC" = "" ] ; then \
		printf "No library found or selected. Maybe PARI/GP is installed in a different directory than what you supplied?\n"; \
		exit; \
	fi; \
	printf "$${LIBLOC%/*} " > ./pari_loc.txt; \
	printf "Searching for \e[32mparipriv.h\e[0m:\n"; \
	FOUNDLOCS=$$(find $$SEARCHLOC -name paripriv.h 2>/dev/null); \
	printf "Here are the places we found \e[32mlibpari.so\e[0m:\n\e[32m"; \
	for poss in $$FOUNDLOCS; do \
		printf "\t$$poss\n"; \
	done; \
	printf "\e[0mYou will now pick which one is correct.\n"; \
	HEADERLOC=""; \
	for poss in $$FOUNDLOCS; do \
		printf "Is \e[35m$$poss\e[0m correct? "; \
		read -r -p "(y/n) " RESPONSE; \
		if [ "$$RESPONSE" = "y" ]; then \
			HEADERLOC=$$poss; \
			break; \
		fi; \
	done; \
	if [ "$$HEADERLOC" = "" ] ; then \
		printf "No library found or selected. Maybe PARI/GP is installed in a different directory than what you supplied?\n"; \
		exit; \
	fi; \
	printf "$${HEADERLOC%/*}" >> ./pari_loc.txt; \
	echo "Setup complete!"
