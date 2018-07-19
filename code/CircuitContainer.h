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

#ifndef CIRCUITCONTAINER_H
#define CIRCUITCONTAINER_H

#include "common.h"
//#include "Circuit.h"
#include "BigIntLib.h"
#include <vector>

class CircuitContainer {
public:
  CircuitContainer();
  ~CircuitContainer();
  void init(uint32 value_size, uint32 random_tape_size, uint32 key_size, uint32 branch_bits, uint32 num_branches, uint32 field_tpye, uint32 party_size);
  void initCipher(uint32 cipher_type);
  void initMiMC();
  void runSign(uchar* x, SignData* sign_data);
  void runVerify(Proof* p, uchar* y, VerifyData* verify_data, uint32 iteration);
  void directEncryption(uchar* x, uchar* y); // Wrapper
  void directMiMC(uchar* x, uchar* y);
  void getParams(uint32* value_size, uint32* random_tape_size, uint32* key_size, uint32* gate_size, uint32* num_view_gates);
  void randomizeKey();
  uchar* getKey();
  uint32 getLastCircuitSignNS();
  uint32 getLastCircuitVerifyNS();
  uint64 getLastDirectCallCycles();
  uint32 getLastDirectCallNS();
  uint32 getCipherNumBranches();
  uint32 getCipherNumMulGates();
  uint32 getCipherNumRounds();
  uint32 getViewSizeBits();
  std::string getCipherName();
  
private:
  std::string cipher_name_;
  uint32 party_size_;
  uint32 value_size_; // Size of x and y (assuming x_size == y_size)
  uint32 random_tape_size_;
  uint32 key_size_; // Size of the key
  uint32 branch_size_; // Part of the gate which is actually used
  uint32 gate_size_; // Size of the gate operations (MUL, ADD, ...)
  uint32 gate_num_words_;
  uint32 branch_bits_;
  uint32 hash_size_;
  uint32 num_rounds_;
  uint32 num_mul_gates_;
  uint32 num_view_gates_;
  uint32 num_feistel_branches_; // Used for Feistel networks
  word** round_constants_; // Used for circuits which require round constants
  uint32 num_round_constants_;
  uchar* value_;
  uchar** value_shares_; // Should always be party_size entries with size value_size
  //word*** value_shares_f_;
  uchar* key_;
  uchar** key_shares_;
  //word*** key_shares_f_;
  word** random_numbers_;
  word* squaring_precomp_;
  std::vector<uint32> feistel_branch_indices_;
  SignData* sign_data_;
  Proof* proof_;
  VerifyData* verify_data_;

  // Circuit inner computation
  uint32 e_; // e chosen be ZKB++ protocol
  uint32 iteration_; // iteration number in the ZKB++ protocol (maybe better to just pass the right data from the outside?)
  uint32 current_mul_gate_;

  // Helpers
  uint32 last_circuit_sign_time_;
  uint32 last_circuit_verify_time_;
  uint64 last_direct_call_cycles_;
  uint32 last_direct_call_time_;

  // Function pointers
  void (CircuitContainer::*direct_function_)(uchar* x, uchar* y);
  void (CircuitContainer::*circuit_function_)(word* value_shares_f, word* key_shares_f, uint32 party_size);
  void (CircuitContainer::*add_c_shared_function_)(word* a_shares, word* b, word* c_shares);
  void (CircuitContainer::*add_shared_function_)(word* a_shares, word* b_shares, word* c_shares);
  void (CircuitContainer::*mul_shared_function_)(word* a_shares, word* b_shares, word* c_shares);
  void (CircuitContainer::*squ_shared_function_)(word* a_shares, word* b_shares, word* c_shares);
  void (CircuitContainer::*squ_shared_experimental_function_)(word* a_shares, word* b_shares, word* c_shares);
  void (CircuitContainer::*cube_shared_function_)(word* a_shares, word* c_shares); // b_shares probably not needed
  void (CircuitContainer::*calc_last_share_function_)(uchar* share_3, uchar* value, uchar* share_1, uchar* share_2); /* DEPRECATED */
  void (CircuitContainer::*fill_value_shares_function_)(word* shares_words, uchar** shares, uint32 party_size); /* DEPRECATED */
  void (CircuitContainer::*prepare_shares_field_sign_function_)(uchar* x, void* key_shares, word* value_shares_f, word* key_shares_f);
  void (CircuitContainer::*prepare_shares_field_verify_function_)(uchar** key_shares, word* value_shares_f, word* key_shares_f);
  void (CircuitContainer::*output_shares_to_bytes_function_)(word* output_shares, uint32 party_size);
  void (CircuitContainer::*verify_calc_last_share_function_)(uchar* y, word* value_shares_f);

  // View allocation (used for SignData* and VerifyData*)
  //View* createView(uint32 num_gates, uint32 gate_size);

  // Specific circuits (all circuits should take the same parameters!)
  void circuitMiMC(word* value_shares_f, word* key_shares_f, uint32 party_size);

  // Circuit preparation and finalization
  void initMatrixM();
  void initMatrixMPrime();
  void beforeSign(uchar* x, SignData* sign_data, word* value_shares_f, word* key_shares_f);
  void afterSign(word* value_shares_f);
  void beforeVerify(Proof* p, uchar* y, VerifyData* verify_data, uchar** key_shares, word* value_shares_f, word* key_shares_f, uint32 iteration);
  void afterVerify(uchar* y, word* value_shares_f);

  // Circuit
  void prepareRoundConstants(uint32 length);
  void destroyRoundConstants();
  void prepareRandomNumbers(uchar** random_tapes, uint32 party_size);
  void destroyRandomNumbers();
  void addSharedSign(word* a_shares, word* b_shares, word* c_shares);
  void addSharedVerify(word* a_shares, word* b_shares, word* c_shares);
  void addCSharedSign(word* a_shares, word* b, word* c_shares);
  void addCSharedVerify(word* a_shares, word* b, word* c_shares);
  void mulSharedSign(word* a_shares, word* b_shares, word* c_shares);
  void mulSharedVerify(word* a_shares, word* b_shares, word* c_shares);
  void squSharedSign(word* a_shares, word* b_shares, word* c_shares);
  void squSharedVerify(word* a_shares, word* b_shares, word* c_shares);
  void squSharedExperimental17Sign(word* a_shares, word* b_shares, word* c_shares);
  void squSharedExperimental17Verify(word* a_shares, word* b_shares, word* c_shares);
  void squSharedExperimental33Sign(word* a_shares, word* b_shares, word* c_shares);
  void squSharedExperimental33Verify(word* a_shares, word* b_shares, word* c_shares);
  void cubeSharedSign(word* a_shares, word* c_shares);
  void cubeSharedVerify(word* a_shares, word* c_shares);
  void copyShares(word* from_shares, word* to_shares, uint32 party_size);

  // Additional calculations
  void calcLastShareDefault(uchar* share_3, uchar* value, uchar* share_1, uchar* share_2); /* DEPRECATED */
  void calcLastShareFeistel(uchar* share_3, uchar* value, uchar* share_1, uchar* share_2); /* DEPRECATED */
  void fillValueSharesDefault(word* shares_words, uchar** shares, uint32 party_size); // /* DEPRECATED */
  void fillValueSharesFeistel(word* shares_words, uchar** shares, uint32 party_size); // /* DEPRECATED */
  void prepareSharesFieldSign(uchar* x, void* key_shares, word* value_shares_f, word* key_shares_f);
  void prepareSharesFieldSign4B(uchar* x, void* key_shares, word* value_shares_f, word* key_shares_f);
  void prepareSharesFieldSignGeneric(uchar* x, void* key_shares, word* value_shares_f, word* key_shares_f);
  void prepareSharesFieldVerify(uchar** key_shares, word* value_shares_f, word* key_shares_f);
  void prepareSharesFieldVerify4B(uchar** key_shares, word* value_shares_f, word* key_shares_f);
  void prepareSharesFieldVerifyGeneric(uchar** key_shares, word* value_shares_f, word* key_shares_f);
  void outputSharesToBytes(word* output_shares, uint32 party_size);
  void outputSharesToBytes4B(word* output_shares, uint32 party_size);
  void outputSharesToBytesGeneric(word* output_shares, uint32 party_size);
  void verifyCalcLastShare(uchar* y, word* value_shares_f);
  void verifyCalcLastShare4B(uchar* y, word* value_shares_f);
  void verifyCalcLastShareGeneric(uchar* y, word* value_shares_f);

  // Util
  uint64 rdtsc();
};

#endif // CIRCUITCONTAINER_H