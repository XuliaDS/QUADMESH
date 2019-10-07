all:
	make -f mkEspQuads.make clean
	#make -f mkNewton.make
	#make -f mkForceQuads.make
	make -f mkServer.make
	make -f mkRegQuads.make
#	make -f mkespQuads.make
#	make -f mkOldQ.make

