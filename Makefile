all:
	make -f mkRegQuads.make clean
	make -f mkForceQuads.make
	make -f mkServer.make
	make -f mkRegQuads.make
