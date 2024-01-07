GDILIB = /usr/lib/w32api/libgdi32.a
LIBRARY_OBJECTS = manipulate_structures7.c

# CC_FLAGS =  -malign-double -fstrict-aliasing -mtune=pentiumpro -mno-cygwin
#CFLAG1 = -mno-cygwin -pg
CFLAG1 = -O2
#CFLAG1= -pg
# CFLAG1 =  -malign-double -O2 -fstrict-aliasing -mtune=pentiumpro -mno-cygwin


main.exe : main.c myres.res manipulate_structures7.o createxyzr.o struct2chain.o new_struct2chain.o jones.o myutils-output.o utils-main.o myvolyme.o onechain.o vetortransfor.o plotstruct.o
	gcc -mwindows $(CFLAG1) -lm  main.c manipulate_structures7.o struct2chain.o createxyzr.o new_struct2chain.o jones.o myutils-output.o utils-main.o myvolyme.o onechain.o vetortransfor.o plotstruct.o myres.res -o $@


myres.res : myres.rc StrchBlt.h 
	windres $< -O coff -o $@

manipulate_structures7.o : manipulate_structures7.c structures.h
	gcc -c $(CFLAG1) manipulate_structures7.c

createxyzr.o : createxyzr.c structures.h
	gcc -c $(CFLAG1) createxyzr.c

struct2chain.o : struct2chain.c structures.h
	gcc -c $(CFLAG1) struct2chain.c

utils-main.o : utils-main.c myutils.h
	gcc -c $(CFLAG1) utils-main.c

myvolyme.o : myvolyme.c myutils.h
	gcc -c $(CFLAG1) myvolyme.c

myutils-output.o : myutils-output.c myutils.h
	gcc -c $(CFLAG1) myutils-output.c

jones.o : jones.c structures.h
	gcc -c $(CFLAG1) jones.c

new_struct2chain.o : new_struct2chain.c structures.h
	gcc -c $(CFLAG1) new_struct2chain.c

onechain.o : onechain.c structures.h
	gcc -c $(CFLAG1) onechain.c

vetortransfor.o : vetortransfor.c myutils.h
	gcc -c $(CFLAG1) vetortransfor.c

plotstruct.o : plotstruct.c structures.h
	gcc -c $(CFLAG1) plotstruct.c