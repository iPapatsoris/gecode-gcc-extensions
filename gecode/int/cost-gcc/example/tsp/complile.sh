g++ -ggdb -I. -fcx-limited-range -fno-signaling-nans -fno-rounding-math -ffinite-math-only -fno-math-errno -fno-strict-aliasing -O3 -fvisibility=hidden -ggdb -std=c++11 -pipe -Wno-unknown-pragmas -Wall -Wextra -fPIC -pthread -DNDEBUG   -DQT_NO_DEBUG -DQT_PRINTSUPPORT_LIB -DQT_WIDGETS_LIB -DQT_GUI_LIB -DQT_CORE_LIB  -I. -isystem /usr/include/x86_64-linux-gnu/qt5 -isystem /usr/include/x86_64-linux-gnu/qt5/QtPrintSupport -isystem /usr/include/x86_64-linux-gnu/qt5/QtWidgets -isystem /usr/include/x86_64-linux-gnu/qt5/QtGui -isystem /usr/include/x86_64-linux-gnu/qt5/QtCore -I. -isystem /usr/include/libdrm -I/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++ \
-c -o cost-gcc-example.o  cost-gcc-example.cpp -pg -g

g++ -ggdb -I. -fcx-limited-range -fno-signaling-nans -fno-rounding-math -ffinite-math-only -fno-math-errno -fno-strict-aliasing -O3 -fvisibility=hidden -ggdb -std=c++11 -pipe -Wno-unknown-pragmas -Wall -Wextra -fPIC -pthread -DNDEBUG   -DQT_NO_DEBUG \
-c -o read-input.o  read-input.cpp -pg -g

g++ -ggdb -I. -fcx-limited-range -fno-signaling-nans -fno-rounding-math -ffinite-math-only -fno-math-errno -fno-strict-aliasing -O3 -fvisibility=hidden -ggdb -std=c++11 -pipe -Wno-unknown-pragmas -Wall -Wextra -fPIC -pthread -DNDEBUG   -DQT_NO_DEBUG -DQT_PRINTSUPPORT_LIB -DQT_WIDGETS_LIB -DQT_GUI_LIB -DQT_CORE_LIB  -I. -isystem /usr/include/x86_64-linux-gnu/qt5 -isystem /usr/include/x86_64-linux-gnu/qt5/QtPrintSupport -isystem /usr/include/x86_64-linux-gnu/qt5/QtWidgets -isystem /usr/include/x86_64-linux-gnu/qt5/QtGui -isystem /usr/include/x86_64-linux-gnu/qt5/QtCore -I. -isystem /usr/include/libdrm -I/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++ \
-c -o flow-graph.o  ../../flow-graph.cpp -pg -g


g++ -ggdb -o cost-gcc cost-gcc-example.o flow-graph.o read-input.o -pg -g -L. -I. -fcx-limited-range -fno-signaling-nans -fno-rounding-math -ffinite-math-only -fno-math-errno -fno-strict-aliasing -O3 -fvisibility=hidden -ggdb -std=c++11 -pipe -Wno-unknown-pragmas -Wall -Wextra -fPIC -pthread -DNDEBUG  \
-lgecodeflatzinc -lgecodedriver -lgecodegist -lgecodesearch -lgecodeminimodel -lgecodeset -lgecodefloat  -lgecodeint -lgecodekernel -lgecodesupport   -lQt5PrintSupport -lQt5Widgets -lQt5Gui -lQt5Core -lGL -lpthread 
