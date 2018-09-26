/*
MIT License

Copyright (c) 2018 Markus Schofnegger

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "common.h"
#include "CircuitContainer.h"

class ZKBPP {
public:
  ZKBPP();
  ~ZKBPP();
  void init(uint32 party_size, uint32 num_iterations, CircuitContainer* cc, bool print_result = true);
  Proof* sign(uchar* x);
  bool verify(Proof* p, uchar* x, uchar* y);
  uint64 getLastGenSignNS();
  uint32 getLastSignNS();
  uint64 getLastGenVerifyNS();
  uint32 getLastVerifyNS();
  uint32 getExpectedSignatureSizeBits();
private:
  uint32 party_size_;
  uint32 num_iterations_;
  CircuitContainer* circuit_;
  ContainerCD* commitment_;

  // Circuit values
  uint32 circuit_value_size_; // Size of x and y (assuming x_size == y_size)
  uint32 circuit_key_size_; // Size of the key
  uint32 circuit_gate_size_; // Size of the gate operations (MUL, ADD, ...)
  uint32 circuit_num_view_gates_; // Number of gates contributing to the View

  // Inner size values
  uint32 hash_size_;
  uint32 random_tape_size_;
  uint32 view_size_;

  // Timing
  uint64 last_gensign_time_;
  uint64 last_sign_time_;
  uint64 last_genverify_time_;
  uint64 last_verify_time_;

  // Utils
  bool print_result_;

  // Allocation and deletion of structs
  View* createView();
  void destroyView(View* v);
  SignData* createSignData(void* sign_views_buffer, void* random_tapes_buffer, uint32 iteration);
  void destroySignData(SignData* sign_data);
  ContainerSignData* createContainerSignData(void* sign_views_buffer, void* random_tapes_buffer);
  void destroyContainerSignData(ContainerSignData* csd);
  ContainerVerifyData* createContainerVerifyData(void* verify_views_buffer);
  void destroyContainerVerifyData(ContainerVerifyData* cvd);
  VerifyData* createVerifyData(void* verify_views_buffer, uint32 iteration);
  void destroyVerifyData(VerifyData* verify_data);
  C* createC();
  void destroyC(C* c);
  D* createD(uint32 k_View_size);
  void destroyD(D* d);
  ContainerCD* createContainerCD(uint32 k_View_size, uint32 party_size); // party_size has to remain (work with -1 for Verify step)
  void destroyContainerCD(ContainerCD* ccd, uint32 party_size); // party_size has to remain (work with -1 for Verify step)
  A* createA();
  void destroyA(A* a);
  ContainerA* createContainerA();
  void destroyContainerA(ContainerA* ca);
  B* createB();
  void destroyB(B* b);
  Z* createZ(bool create_x_3);
  void destroyZ(Z* z);
  Proof* createProof();
  void destroyProof(Proof* p);

  // Struct content
  void fillCDSign(ContainerCD* ccd, uint32 iteration, SignData* sign_data, uchar* hash_data);
  void fillCDVerify(ContainerCD* ccd, uint32 iteration, Proof* p, VerifyData* verify_data, uchar* hash_data);
  void fillASign(ContainerA* ca, uint32 iteration, SignData* sign_data, ContainerCD* ccd);
  void fillAVerify(ContainerA* ca, uint32 iteration, Proof* p, VerifyData* verify_data, ContainerCD* ccd);
  void fillProof(Proof* p, uint32 iteration, SignData* sign_data, ContainerCD* ccd);

  // Signature
  void buildChallengeHash(uchar* destination, ContainerA* ca);

  // Hashs
  void SHA256(uchar* destination, uchar* data, uint32 data_size);
  void SHA256Prime(uchar* destination, uchar* data, uint32 data_size);
  void SHA256Dash(uchar* destination, uchar* data, uint32 data_size);

  // Other
  void printView(uchar* view);
  void printProof(Proof* p, bool print_view);
  void printDataAsHex(uchar* data, uint32 data_size, bool format);
};
