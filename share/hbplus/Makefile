##############################################################################
#Choose your C compiler.  Usually gcc.
#CC     = cc
CC    = gcc

#Also the flags . . .
# CFLAGS = -g -DBSM -fullwarn
CFLAGS = -O3

#Sometimes, you need to include this line if HBPLUS does not compile.
#COBJECT = /usr/lib/libm.a
COBJECT = -lm

#SHELL = /usr/bin/csh
##############################################################################

##############################################################################
#To ensure you compile all four . . 
programs :hbplus
##############################################################################




##############################################################################
hbplus : hbp_gen.o hbp_inpdb.o hbp_findh.o hbp_hhb.o hbp_qnh.o hbp_main.o ;
	$(CC) $(CFLAGS) -o hbplus \
	hbp_gen.o hbp_inpdb.o hbp_findh.o hbp_hhb.o hbp_qnh.o hbp_main.o $(COBJECT)

hbp_gen.o   : hbplus.h
	$(CC) -c $(CFLAGS) hbp_gen.c
hbp_inpdb.o : hbplus.h
	$(CC) -c $(CFLAGS) hbp_inpdb.c
hbp_findh.o : hbplus.h
	$(CC) -c $(CFLAGS) hbp_findh.c
hbp_hhb.o   : hbplus.h
	$(CC) -c $(CFLAGS) hbp_hhb.c 
hbp_dom.o   : hbplus.h
	$(CC) -c $(CFLAGS) hbp_dom.c 
hbp_qnh.o   : hbplus.h
	$(CC) -c $(CFLAGS) hbp_qnh.c 
hbp_main.o  : hbplus.h
	$(CC) -c $(CFLAGS) hbp_main.c
##############################################################################


