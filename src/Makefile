CC = icpc -O3
CFLAGS = -I${MKLROOT}/include
LD_FLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl


target = dqmc
src = $(wildcard *.cpp)
obj = $(patsubst %.cpp, %.o, $(src))

$(target) : $(obj)
	$(CC) $(LD_FLAGS) $(obj) -o $(target)

main.o: main.cpp RingHubbard.h SquareHubbard.h cfg.h dqmc.h model.h
	$(CC) -c $<

dqmc.o: dqmc.cpp dqmc.h model.h gfun.h gtau.h measure.h
	$(CC) $(CFLAGS) -c $<

cfg.o: cfg.cpp cfg.h
	$(CC) -c $<

gfun.o: gfun.cpp gfun.h matb.h model.h seqb.h
	$(CC) $(CFLAGS) -c $<

gtau.o: gtau.cpp gtau.h matb.h model.h
	$(CC) $(CFLAGS) -c $<

matb.o: matb.cpp matb.h model.h
	$(CC) $(CFLAGS) -c $<

measure.o: measure.cpp measure.h gtau.h
	$(CC) $(CFLAGS) -c $<

RingHubbard.o: RingHubbard.cpp RingHubbard.h model.h
	$(CC) -c $<

SquareHubbard.o: SquareHubbard.cpp SquareHubbard.h model.h
	$(CC) -c $<

seqb.o: seqb.cpp seqb.h matb.h
	$(CC) $(CFLAGS) -c $<

util.o: util.cpp util.h
	$(CC) $(CFLAGS) -c $<


.PHONY: clean
clean:
	rm -rf $(obj) $(target)