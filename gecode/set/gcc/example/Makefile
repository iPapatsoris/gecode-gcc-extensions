OBJS =  sym-gcc-example.o read-input.o flow-graph.o

CC 	= g++
FLAGS	= -fcx-limited-range -fno-signaling-nans -fno-rounding-math -ffinite-math-only -fno-math-errno -fno-strict-aliasing -O3 -fvisibility=hidden -std=c++11 -pipe -Wno-unknown-pragmas -Wall -Wextra -fPIC -pthread -DNDEBUG
SRCFLAGS	= $(FLAGS) -DQT_NO_DEBUG 
LIBS	= -lgecodedriver -lgecodegist -lgecodesearch -lgecodeminimodel -lgecodeint -lgecodeset -lgecodekernel -lgecodesupport -lpthread 

all: sym-gcc gendata

sym-gcc: $(OBJS) $(HEADER)
	$(CC) -o sym-gcc $(OBJS) -L. -I. $(FLAGS) $(LIBS)

sym-gcc-example.o: sym-gcc-example.cpp
	$(CC) $(SRCFLAGS) -I. -isystem /usr/include/libdrm -c -o sym-gcc-example.o sym-gcc-example.cpp

read-input.o: read-input.cpp
	$(CC) $(SRCFLAGS) -c -o read-input.o read-input.cpp

flow-graph.o: ../flow-graph.cpp
	$(CC) $(SRCFLAGS) -c -o flow-graph.o ../flow-graph.cpp

gendata: gendata.cpp
	g++ gendata.cpp -o gendata

clean:
	rm -f sym-gcc gendata $(OBJS)

count:
	wc -l *.cpp *.hpp ../*.cpp ../*.hpp 