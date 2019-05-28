CFLAGS+="-I/usr/local/Cellar/openssl/1.0.2q/include -L/usr/local/Cellar/openssl/1.0.2q/lib/"
${CXX} -std=c++0x -g -O2 -march=native -mtune=native ${CFLAGS} main.cpp affinity_osx.cpp ZKBPP.cpp CircuitContainer.cpp BigIntLib.cpp -o zkbpp_test -lcrypto -lm
