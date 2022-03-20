CC = g++ -O3
CFLAGS = -m64 -I/opt/intel/mkl/include
LD_FLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

all : modelprint matexp checkmatb gstab checkgtau checkcfg LUfact

modelprint: modelprint.o RingHubbard.o SquareHubbard.o
	$(CC) $^ -o $@

modelprint.o: modelprint.cpp ../RingHubbard.h ../SquareHubbard.h
	$(CC) -c $<

RingHubbard.o: ../RingHubbard.cpp ../RingHubbard.h
	$(CC) -c $<

SquareHubbard.o : ../SquareHubbard.cpp ../SquareHubbard.h
	$(CC) -c $<

#==========================================================

matexp : matexp.o util.o
	$(CC) $(LD_FLAGS) $^ -o $@

matexp.o: matexp.cpp ../util.h
	$(CC) $(CFLAGS) -c $<

util.o: ../util.cpp ../util.h
	$(CC) $(CFLAGS) -c $<

#==========================================================

checkmatb: checkmatb.o RingHubbard.o matb.o util.o
	$(CC) $(LD_FLAGS) $^ -o $@

checkmatb.o: checkmatb.cpp ../RingHubbard.h ../matb.h ../util.h
	$(CC) $(CFLAGS) -c $<

matb.o: ../matb.cpp ../matb.h ../util.h
	$(CC) $(CFLAGS) -c $<

#==========================================================

gstab: gstab.o util.o seqb.o matb.o RingHubbard.o gfun.o
	$(CC) $(LD_FLAGS) $^ -o $@

gstab.o: gstab.cpp ../gfun.h ../seqb.h ../matb.h ../util.h ../RingHubbard.h
	$(CC) $(CFLAGS) -c $<

seqb.o: ../seqb.cpp ../seqb.h ../util.h
	$(CC) $(CFLAGS) -c $<

gfun.o: ../gfun.cpp ../gfun.h ../util.h
	$(CC) $(CFLAGS) -c $<

#==========================================================

checkgtau: checkgtau.o gtau.o gfun.o util.o RingHubbard.o matb.o seqb.o
	$(CC) $(LD_FLAGS) $^ -o $@

checkgtau.o: checkgtau.cpp ../gtau.h ../gfun.h ../util.h ../RingHubbard.h
	$(CC) $(CFLAGS) -c $<

gtau.o: ../gtau.cpp ../gtau.h ../util.h ../matb.h
	$(CC) $(CFLAGS) -c $<

#==========================================================

checkcfg: checkcfg.o cfg.o
	$(CC) $^ -o $@

checkcfg.o: checkcfg.cpp ../cfg.h
	$(CC) -c $<

cfg.o: ../cfg.cpp ../cfg.h
	$(CC) -c $<

#==========================================================

LUfact: LUfact.o util.o
	$(CC) $(LD_FLAGS) $^ -o $@

LUfact.o: LUfact.cpp
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf *.o modelprint matexp checkmatb gstab checkgtau checkcfg LUfact