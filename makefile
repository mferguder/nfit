_OBJS =	boxchain.o dataset.o funclmdif.o funsupport.o globalVariables.o \
        modelcalculator.o nrutil.o Para.o toad.o toadcmd.o toadmisc.o \
        tvds.o tvDSfit.o tvImg.o tvLinfitDriver.o utable.o \
        fileTools.o interp2d.o interp2d_spline.o bicubic.o \
        nfit.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

CC = g++
ODIR = obj
SDIR = src
IDIR = inc

TCL = tcl8.4
TK = tk8.4
BLT = BLT24

CFLAGS = -Wall -O3 -ffast-math -fPIC -lpthread
CLINK = -lm -l$(TCL) -l$(TK) -ltiff -l$(BLT) -lpthread -lmydll \
        -lgsl -lgslcblas -Wl,-rpath,.
CPATH = -I/usr/include/$(TCL) -I./$(IDIR) -L. -I.

all: libtoad toad

toad: libtoad
	gcc $(ODIR)/toad.o $(CPATH) -ltoad_threaded $(CLINK) -o toad 

libtoad: $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) $(CPATH) $(CLINK) -shared -o libtoad_threaded.so 

$(ODIR)/%.o: $(SDIR)/%.cpp $(IDIR)/%.h
	$(CC) -c $(CPATH) $(CFLAGS) -o $@ $< 

clean:
	rm -f $(ODIR)/*.o libtoad_threaded.so


	

