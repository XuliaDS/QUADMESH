all:
	make -f mkXQuads.make clean
	#make -f mkForceQuads.make
	make -f mkServer.make
	make -f mkXQuads.make
#	make -f mkespQuads.make
#	make -f mkOldQ.make

