OBJS =  cost-gcc-example.o read-input.o flow-graph.o

CC 	= g++
FLAGS	= -fcx-limited-range -fno-signaling-nans -fno-rounding-math -ffinite-math-only -fno-math-errno -fno-strict-aliasing -O3 -fvisibility=hidden -std=c++11 -pipe -Wno-unknown-pragmas -Wall -Wextra -fPIC -pthread -DNDEBUG
SRCFLAGS	= $(FLAGS) -DQT_NO_DEBUG 
LIBS	= -lgecodedriver -lgecodegist -lgecodesearch -lgecodeminimodel -lgecodeint -lgecodekernel -lgecodesupport -lpthread 

all: $(OBJS) $(HEADER)
	$(CC) -o cost-gcc $(OBJS) -L. -I. $(FLAGS) $(LIBS)

cost-gcc-example.o: cost-gcc-example.cpp
	$(CC) $(SRCFLAGS) -I. -isystem /usr/include/libdrm -c -o cost-gcc-example.o cost-gcc-example.cpp

read-input.o: read-input.cpp
	$(CC) $(SRCFLAGS) -c -o read-input.o read-input.cpp

flow-graph.o: ../flow-graph.cpp
	$(CC) $(SRCFLAGS) -c -o flow-graph.o ../flow-graph.cpp

clean:
	rm -f cost-gcc $(OBJS)

count:
	wc -l *.cpp *.hpp ../*.cpp ../*.hpp 