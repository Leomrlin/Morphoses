
# INTV = -DLONG
INTE = -DEDGELONG

#compilers
$(info ************  Using CILK ************)
PCC = g++-7
PCFLAGS = -std=c++14 -fcilkplus -lcilkrts -g -O3 -DCILK $(INTV) $(INTE)

# dependencies
COMMON = ../../core/common/utils.h ../../core/common/parseCommandLine.h ../../core/common/parallel.h ../../core/common/quickSort.h ../../core/common/blockRadixSort.h ../../core/common/transpose.h ../../core/common/rwlock.h
LOCAL_COMMON = ../common/graphIO.h ../common/fileUtils.h

CONVERTERS = SNAPtoAdjConverter AddWeightToSNAPFile

.PHONY: all clean

all: $(CONVERTERS)

% : %.C $(COMMON) $(LOCAL_COMMON)
	$(PCC) $(PCFLAGS) -o $@ $<

clean :
	rm -f *.o $(CONVERTERS)
