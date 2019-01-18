#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
#ifdef ESP_BLOC
ODIR = obj
TDIR = BUILDS
#else
#ODIR = .
#TDIR = $(ESP_ROOT)/bin
#endif
FILE = quadServer
SFILE = $(FILE).c
SDIR = SRC
$(TDIR)/$(FILE):	$(ODIR)/$(FILE).o  $(LDIR)/libwsserver.a
	$(CC) -g -o $(TDIR)/$(FILE) $(ODIR)/$(FILE).o \
		-L$(LDIR) -lwsserver -legads -lpthread -L/usr/local/lib -lnlopt -lz $(RPATH) -lm

$(ODIR)/$(FILE).o:	$(SDIR)/$(SFILE) $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -g -c $(COPTS) $(DEFINE) -I$(IDIR) $(SDIR)/$(SFILE) -o $(ODIR)/$(FILE).o

clean:
	-rm $(ODIR)/$(FILE).o

cleanall:	clean
	-rm $(TDIR)/$(FILE)
