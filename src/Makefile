CC = mpic++


# -ansi to be compatible with ansi C
# -Wall: display all warning message
# -g: add debug information to object file. to enable gdb debugging
# -I: add including directory
# -O6: optimization options
#  -DNDEBUG: disable asserts
# 
CCFLAGS = -ansi -Wall -O6 -DNDEBUG -I../../lib

#CCFLAGS = -ansi -Wall -O6 -g -I../../lib
#CCFLAGS = -ansi -Wall -O6 -g -DALLSTATS
LDFLAGS = 
LOADLIBES = -lm        # the name used by the automatic linker rule
MLBFLAGS = -DMLB

SRCS = main.cpp parser_gr.cc timer.cc parser_ss.cc request.cpp
HDRS = nodearc.h
CODES = sp.exe spC.exe

all: $(CODES)

sp.exe: $(SRCS) $(HDRS)
	$(CC) $(CCFLAGS) -o sp.exe $(SRCS) $(LOADLIBES)

# sqC and mbpC calculate checksum after solving problem.
spC.exe: $(SRCS) $(HDRS)
	$(CC) $(CCFLAGS) -DCHECKSUM -o spC.exe $(SRCS) $(LOADLIBES)

clean:
	rm -f *~ sq.exe mbp.exe sqC.exe mbpC.exe
