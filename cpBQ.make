#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
ODIR  = obj
TDIR  = BUILDS
FILE  = regQuadsB
SFILE = $(FILE).c
SDIR  = SRC
DBG   = #-g -pg
$(TDIR)/prevQ:	$(ODIR)/prevQ.o 
	$(CC) $(DBG) -o $(TDIR)/prevQ $(ODIR)/prevQ.o \
	-L$(LDIR) -legads -lm

$(ODIR)/prevQ.o:	$(SDIR)/$(SFILE) $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h $(SDIR)/regQuads.h
	$(CC) $(DBG)  -c $(COPTS) $(DEFINE) -DSTANDALONE -I$(IDIR) \
	$(SDIR)/$(SFILE) -o $(ODIR)/prevQ.o

clean:
	-rm $(ODIR)/$prevQ.o

cleanall:	clean
	-rm $(TDIR)/prevQ
