# g++ -std=c++0x -g -O2 main.cpp ZKBPP.cpp CircuitContainerCustomMath.cpp BigIntLib.cpp -o zkbpp_test -lcrypto -lm
ulimit -s 32768 # increase stack size for large view sizes (e.g. MiMC)
echo "--- MiMC-256 (Prime field) ---"
valgrind --max-stackframe=33554432 --main-stacksize=33554432 ./zkbpp_test 256 1 0 1 0
echo "--- MiMC-272 (Prime field) ---"
valgrind --max-stackframe=33554432 --main-stacksize=33554432 ./zkbpp_test 272 1 0 1 0