CC = gcc
LD = gcc
CFLAGS = -O3 -Wall -funroll-loops
RM = /bin/rm -f
OBJS = galsim.o graphics/graphics.o 
EXECUTABLE = galsim
LDFLAGS=-L/opt/X11/lib -lX11 -lm

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(EXECUTABLE)
	
graphics/graphics.o: graphics/graphics.c graphics/graphics.h
	make -C graphics

galsim.o: galsim.c
	$(CC) $(CFLAGS) -c galsim.c

clean:
	$(RM) $(EXECUTABLE) $(OBJS) && make -C graphics clean