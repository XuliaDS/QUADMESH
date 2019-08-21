#
IDIR  = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR  = $(ESP_ROOT)/lib
ODIR  = obj
TDIR  = BUILDS
FILE  = minL
SFILE = $(FILE).cpp
SDIR  = SRC
DBG   = -g -pg
$(TDIR)/$(FILE):	$(ODIR)/$(FILE).o 
	$(CXX) $(DBG) -o $(TDIR)/$(FILE) $(ODIR)/$(FILE).o \
	-L$(LDIR) -legads -lm

$(ODIR)/$(FILE).o:	$(SDIR)/$(SFILE) $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h
	$(CXX) -c $(CPPOPT) $(DEFINE) -I$(IDIR) -I../include \
	$< -o $@

clean:
	-rm $(ODIR)/$(FILE).o

cleanall:	clean
	-rm $(TDIR)/$(FILE)
