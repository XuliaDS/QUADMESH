#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
ODIR  = obj
TDIR  = BUILDS
FILE  = espregq
SFILE = $(FILE).c
SDIR  = SRC
$(TDIR)/$(FILE):	$(ODIR)/$(FILE).o 
	$(CC) -g -pg -o $(TDIR)/$(FILE) $(ODIR)/$(FILE).o \
		-L$(LDIR) -legads -lm

$(ODIR)/$(FILE).o:	$(SDIR)/$(SFILE) $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(SDIR)/$(FILE).h
	$(CC) -g -pg -c $(COPTS) $(DEFINE) -DSTANDALONE -I$(IDIR) \
		 $(SDIR)/$(SFILE) -o $(ODIR)/$(FILE).o

clean:
	-rm $(ODIR)/$(FILE).o

cleanall:	clean
	-rm $(TDIR)/$(FILE)
