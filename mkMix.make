#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
ODIR  = obj
TDIR  = BUILDS
FILE  = mixq
SFILE = $(FILE).c
SDIR  = SRC
DBG   = -g -pg
$(TDIR)/$(FILE):	$(ODIR)/$(FILE).o 
	$(CC) $(DBG) -o $(TDIR)/$(FILE) $(ODIR)/$(FILE).o \
	-L$(LDIR) -legads -lm

$(ODIR)/$(FILE).o:	$(SDIR)/$(SFILE) $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(SDIR)/regQuads.h
	$(CC) $(DBG)  -c $(COPTS) $(DEFINE) -DSTANDALONE -I$(IDIR) \
	$(SDIR)/$(SFILE) -o $(ODIR)/$(FILE).o

clean:
	-rm $(ODIR)/$(FILE).o

cleanall:	clean
	-rm $(TDIR)/$(FILE)
