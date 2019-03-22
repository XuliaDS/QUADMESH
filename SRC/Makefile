#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
#ifdef ESP_BLOC
ODIR = .
TDIR = .
#else
#ODIR = .
#TDIR = $(ESP_ROOT)/bin
#endif
FILE = regQuads
SFILE = $(FILE).c
SDIR = .

$(TDIR)/$(FILE):	$(ODIR)/$(FILE).o  $(LDIR)/libwsserver.a
	$(CC) -g -o $(TDIR)/$(FILE) $(ODIR)/$(FILE).o \
		-L$(LDIR) -legads $(RPATH) -lm

$(ODIR)/$(FILE).o:	$(SDIR)/$(SFILE) $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(IDIR)/wsserver.h
	$(CC) -g -c $(COPTS) $(DEFINE) -DSTANDALONE -I$(IDIR) \
		$(SDIR)/$(SFILE) -o $(ODIR)/$(FILE).o

clean:
	-rm $(ODIR)/$(FILE).o

cleanall:	clean
	-rm $(TDIR)/$(FILE)
