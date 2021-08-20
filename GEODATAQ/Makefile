# 
# 	MAKEFILE FOR GEODATAQ
#
#

CC = gcc -g
#
# OBJECT MODULES 
#
OBJECTS = geodataqmain.o geodataqinit.o geodataqquery.o geominterp.o Stdio.o

.PHONY:	all tags clean cleanall

TARGET	=	geodataq


geodata:	$(OBJECTS)
	$(CC) $(OBJECTS) -lm -o geodataq 

geodatamain.o: geodatamain.c 
geodataqinit.o: geodataqinit.c Stdio.c geodataq.h Stdio.h
Stdio.o: Stdio.c Stdio.h
geodataqquery.o: geodataqquery.c Stdio.c geodataq.h geominterp.h Stdio.h
geominterp.o: geominterp.c geominterp.h geodataq.h

all: $(TARGET)

clean:
	rm -rf $(OBJECTS)  core.* *.o *~ 

TAGS:
	etags *.c *.h

tags:	TAGS

