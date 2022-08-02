g++ -I. -fcx-limited-range -fno-signaling-nans -fno-rounding-math -ffinite-math-only -fno-math-errno -fno-strict-aliasing -O3 -fvisibility=hidden -std=c++11 -pipe -Wno-unknown-pragmas -Wall -Wextra -fPIC -pthread -DNDEBUG -DQT_NO_DEBUG -I. -isystem /usr/include/libdrm \
-c -o cost-gcc-example.o  cost-gcc-example.cpp -g

g++ -I. -fcx-limited-range -fno-signaling-nans -fno-rounding-math -ffinite-math-only -fno-math-errno -fno-strict-aliasing -O3 -fvisibility=hidden -std=c++11 -pipe -Wno-unknown-pragmas -Wall -Wextra -fPIC -pthread -DNDEBUG   -DQT_NO_DEBUG \
-c -o read-input.o  read-input.cpp -g

g++ -I. -fcx-limited-range -fno-signaling-nans -fno-rounding-math -ffinite-math-only -fno-math-errno -fno-strict-aliasing -O3 -fvisibility=hidden -std=c++11 -pipe -Wno-unknown-pragmas -Wall -Wextra -fPIC -pthread -DNDEBUG   -DQT_NO_DEBUG \
-c -o flow-graph.o  ../flow-graph.cpp -g

g++ -o cost-gcc cost-gcc-example.o flow-graph.o read-input.o -g -L. -I. -fcx-limited-range -fno-signaling-nans -fno-rounding-math -ffinite-math-only -fno-math-errno -fno-strict-aliasing -O3 -fvisibility=hidden -std=c++11 -pipe -Wno-unknown-pragmas -Wall -Wextra -fPIC -pthread -DNDEBUG  \
-lgecodedriver -lgecodegist -lgecodesearch -lgecodeminimodel -lgecodeint -lgecodekernel -lgecodesupport -lpthread 
