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

// Shared implementation of circuits.
#include "CircuitContainer.h"
#include <openssl/rand.h> // For RAND_bytes(.)
#include <openssl/sha.h> // FOR SHA256_DIGEST_LENGTH
#include <cmath> // log2(.)
#include <cstring> // memcpy
#include <vector>
#include <algorithm>
#include <stdint.h> // rdtsc

//#define VERBOSE

CircuitContainer::CircuitContainer() {

}

CircuitContainer::~CircuitContainer() {
  if(this->value_shares_ != NULL) {
    for(uint32 i = 0; i < this->party_size_; i++) {
      delete[] this->value_shares_[i];
    }
    delete[] this->value_shares_;
    this->value_shares_ = NULL;
  }
  if(this->value_ != NULL) delete[] this->value_;
  if(this->key_shares_ != NULL) {
    for(uint32 i = 0; i < this->party_size_; i++) {
      delete[] this->key_shares_[i];
    }
    delete[] this->key_shares_;
    this->key_shares_ = NULL;
  }
  if(this->key_ != NULL) delete[] this->key_;

  if(this->round_constants_ != NULL) {
    for(uint32 i = 0; i < this->num_round_constants_; i++) {
      delete[] this->round_constants_[i];
    }
    delete[] this->round_constants_;
  }

  if(this->intermediate_results_ != NULL) {
    for(uint32 i = 0; i < this->num_intermediate_results_; i++) {
      delete[] this->intermediate_results_[i];
    }
    delete[] this->intermediate_results_;
  }

  // Clean up random numbers
  this->destroyRandomNumbers();

  // Clean up matrix M memory (if num_feistel_branches_ > 1)
  
  if(this->num_feistel_branches_ > 1) {
    if(this->matrix_1_ != NULL) delete[] this->matrix_1_;
    if(this->matrix_2_ != NULL) delete[] this->matrix_2_;
  }

  // Clean up numbers for faster squaring
  if(this->squaring_precomp_ != NULL) delete[] this->squaring_precomp_;

  // Clean up BigIntLib
  BigIntLib::CleanUp();
}

void CircuitContainer::init(uint32 value_size, uint32 random_tape_size, uint32 key_size, uint32 branch_bits, uint32 num_branches, uint32 field_type, uint32 party_size) {
  this->cipher_name_ = "";
  this->value_size_ = value_size;
  this->random_tape_size_ = random_tape_size;
  this->key_size_ = key_size;
  this->branch_size_ = ceil(((float)branch_bits / WORD_SIZE) * 8);
  this->gate_size_ = ceil((float)branch_bits / 64) * 8;
  this->gate_num_words_ = ceil(float(this->gate_size_) / (WORD_SIZE / 8));
  this->branch_bits_ = branch_bits;
  this->hash_size_ = SHA256_DIGEST_LENGTH;
  this->party_size_ = party_size;
  this->value_ = new uchar[value_size];
  memset(this->value_, 0, value_size);
  this->value_shares_ = new uchar*[this->party_size_];
  this->key_ = new uchar[key_size];
  this->key_shares_ = new uchar*[this->party_size_];
  //RAND_bytes(this->key_, key_size);
  this->randomizeKey();
  for(uint32 i = 0; i < party_size; i++) {
    this->value_shares_[i] = new uchar[value_size];
    memset(this->value_shares_[i], 0, value_size);
    this->key_shares_[i] = new uchar[key_size];
    memset(this->key_shares_[i], 0, key_size);
  }
  this->matrix_1_ = NULL;
  this->matrix_2_ = NULL;

  this->round_constants_ = NULL;
  this->intermediate_results_ = NULL;
  this->squaring_precomp_ = NULL;
  this->num_feistel_branches_ = num_branches;
  this->feistel_branch_indices_.push_back(0); // Push back single index for ciphers without feistel branches

  // Init BigIntLib
  BigIntLib::Init(branch_bits, field_type);
}

void CircuitContainer::initCipher(uint32 cipher_type) {
  if(cipher_type == 1) {
    this->cipher_name_ = "MiMC";
    this->initMiMC();
  }
  else {
    std::cout << "Error: Unknown cipher type '" << cipher_type << "'." << std::endl;
    exit(1);
  }
}

void CircuitContainer::initMiMC() {
  // Settings
  this->num_rounds_ = ceil((float)(BigIntLib::field_size_bits_) / log2(3));
  this->num_mul_gates_ = 2 + (this->num_rounds_ * 2); // In MiMC x^3 at the beginning and then one x^3 for each round; every x^3 consists of two MUL gates
  this->num_intermediate_results_ = this->num_mul_gates_;
  //this->num_view_gates_ = this->num_mul_gates_ + 1; // + 1 for last gate (y share)
  this->num_view_gates_ = this->num_mul_gates_;
  this->num_round_constants_ = this->num_rounds_;

  // Round constants
  this->prepareRoundConstants(this->num_round_constants_);

  // Function pointers
  this->direct_function_ = &CircuitContainer::directMiMC;
  this->circuit_function_ = &CircuitContainer::circuitMiMC;
  this->prepare_shares_field_sign_function_ = &CircuitContainer::prepareSharesFieldSign;
  this->prepare_shares_field_verify_function_ = &CircuitContainer::prepareSharesFieldVerify;
  this->output_shares_to_bytes_function_ = &CircuitContainer::outputSharesToBytes;
  this->verify_calc_last_share_function_ = &CircuitContainer::verifyCalcLastShare;

  // Allocate space for random numbers
  this->random_numbers_ = new word*[this->party_size_];
  for(uint32 i = 0; i < this->party_size_; i++) {
    this->random_numbers_[i] = new word[this->num_mul_gates_ * this->gate_num_words_];
    memset(this->random_numbers_[i], 0, this->num_mul_gates_ * this->gate_num_words_ * (WORD_SIZE / 8));
  }

  // Allocate space for intermediate results
  this->intermediate_results_ = new word*[this->num_intermediate_results_];
  for(uint32 i = 0; i < this->num_intermediate_results_; i++) {
    this->intermediate_results_[i] = new word[this->gate_num_words_];
    memset(this->intermediate_results_[i], 0x0, this->gate_size_);
  }

}

void CircuitContainer::runSign(uchar* x, SignData* sign_data) {
  // Stack storage for field values
  //word* test = new word[this->num_feistel_branches_ * this->party_size_ * this->gate_num_words_ * 2];
  alignas(8) word value_shares_f[this->num_feistel_branches_][this->party_size_][this->gate_num_words_];
  alignas(8) word key_shares_f[this->num_feistel_branches_][this->party_size_][this->gate_num_words_];
  memset(value_shares_f, 0, this->num_feistel_branches_ * this->party_size_ * this->gate_size_);
  memset(key_shares_f, 0, this->num_feistel_branches_ * this->party_size_ * this->gate_size_);
  
  // Prepare (should be the same routine for each circuit!)
  //auto circuit_sign_start = std::chrono::high_resolution_clock::now();
  this->beforeSign(x, sign_data, (word*)value_shares_f, (word*)key_shares_f);
  //auto circuit_sign_stop = std::chrono::high_resolution_clock::now();

  // Run circuit
  auto circuit_sign_start = std::chrono::high_resolution_clock::now();
  (this->*circuit_function_)((word*)value_shares_f, (word*)key_shares_f, this->party_size_);
  auto circuit_sign_stop = std::chrono::high_resolution_clock::now();
  //this->last_circuit_sign_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(circuit_sign_stop - circuit_sign_start).count();

  // Clean up (should be the same routine for each circuit!)
  //auto circuit_sign_start = std::chrono::high_resolution_clock::now();
  this->afterSign((word*)value_shares_f);
  //auto circuit_sign_stop = std::chrono::high_resolution_clock::now();

  this->last_circuit_sign_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(circuit_sign_stop - circuit_sign_start).count();
}

void CircuitContainer::runVerify(Proof* p, uchar* x, uchar* y, VerifyData* verify_data, uint32 iteration) {
  // Prepare (should be the same routine for each circuit!)
  // Stack storage for field values
  alignas(8) word value_shares_f[this->num_feistel_branches_][this->party_size_ - 1][this->gate_num_words_];
  alignas(8) word key_shares_f[this->num_feistel_branches_][this->party_size_ - 1][this->gate_num_words_];
  memset(value_shares_f, 0, this->num_feistel_branches_ * (this->party_size_ - 1) * this->gate_size_);
  memset(key_shares_f, 0, this->num_feistel_branches_ * (this->party_size_ - 1) * this->gate_size_);

  #ifdef VERBOSE
  std::cout << "[VERIFY] Entering Verify step with y: ";
  BigIntLib::Print(y, this->value_size_);
  #endif

  uchar* key_shares[this->party_size_ - 1];
  this->beforeVerify(p, x, y, verify_data, key_shares, (word*)value_shares_f, (word*)key_shares_f, iteration);

  // Run circuit
  auto circuit_verify_start = std::chrono::high_resolution_clock::now();
  (this->*circuit_function_)((word*)value_shares_f, (word*)key_shares_f, this->party_size_ - 1);
  auto circuit_verify_stop = std::chrono::high_resolution_clock::now();
  this->last_circuit_verify_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(circuit_verify_stop - circuit_verify_start).count();

  // Clean up (should be the same routine for each circuit!)
  this->afterVerify(y, (word*)value_shares_f);
}

void CircuitContainer::directEncryption(uchar* x, uchar* y) {
  auto direct_call_start = std::chrono::high_resolution_clock::now();
  (this->*direct_function_)(x, y);
  auto direct_call_stop = std::chrono::high_resolution_clock::now();
  this->last_direct_call_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(direct_call_stop - direct_call_start).count();
}

void CircuitContainer::directMiMC(uchar* x, uchar* y) {

  this->current_intermediate_result_ = 0;

  uint64 cycles_begin = this->rdtsc();

  // Key
  word key[this->gate_num_words_];
  memset(key, 0, (this->gate_num_words_) * (WORD_SIZE / 8));
  memcpy(key, this->key_, this->key_size_);

  word temp[this->gate_num_words_];
  memset(temp, 0, (this->gate_num_words_) * (WORD_SIZE / 8));

  // Input and output
  word input[this->gate_num_words_];
  memset(input, 0, (this->gate_num_words_) * (WORD_SIZE / 8));
  memcpy(input, x, this->value_size_);

  // Keep in field
  BigIntLib::TryReduce(input);
  BigIntLib::TryReduce(key);

  word output[this->gate_num_words_];
  memset(output, 0, (this->gate_num_words_) * (WORD_SIZE / 8));

  // ADD with round key (omit adding constant, since c_0 = 0)
  BigIntLib::Add(output, input, key);

  // x^3 calculation
  BigIntLib::Mul(temp, output, output);
  BigIntLib::Mul(input, temp, output);

  // --- New ZKB++ optimization
  memcpy(this->intermediate_results_[this->current_intermediate_result_++], temp, this->gate_size_);
  memcpy(this->intermediate_results_[this->current_intermediate_result_++], input, this->gate_size_);
  // ---

  for(uint32 i = 1; i < this->num_rounds_; i++) {
    // ADD with round key
    BigIntLib::Add(output, input, key);
  
    // ADD with round constant
    BigIntLib::Add(output, output, this->round_constants_[i]);
  
    // x^3 calculation
    BigIntLib::Mul(temp, output, output);
    BigIntLib::Mul(input, temp, output);
    // --- New ZKB++ optimization
    memcpy(this->intermediate_results_[this->current_intermediate_result_++], temp, this->gate_size_);
    memcpy(this->intermediate_results_[this->current_intermediate_result_++], input, this->gate_size_);
    // ---
  }
  // ADD with round key
  BigIntLib::Add(output, input, key);
  //std::cout << "[DIRECT] Value: " << BigIntLib::ToString(output) << std::endl;

  this->last_direct_call_cycles_ = this->rdtsc() - cycles_begin;

  std::cout << "Direct output: " << BigIntLib::ToString(output) << std::endl;
  memcpy(y, output, this->value_size_);
}

void CircuitContainer::getParams(uint32* value_size, uint32* random_tape_size, uint32* key_size, uint32* gate_size, uint32* num_view_gates) {
  //uint32 num_mul_gates = 2 + (this->num_rounds_ * 2); // In MiMC x^3 at the beginning and then one x^3 for each round; every x^2 consists of two MUL gates
  //uint32 num_view_gates = num_mul_gates + 1; // + 1 for last gate (y share)
  //uint32 num_view_gates = num_mul_gates;
  //return num_view_gates;
  *value_size = this->value_size_;
  *random_tape_size = this->random_tape_size_;
  *key_size = this->key_size_;
  *gate_size = this->gate_size_;
  *num_view_gates = this->num_view_gates_;
}

void CircuitContainer::randomizeKey() {
  RAND_bytes(this->key_, this->key_size_);
  //memset(this->key_, 0xAB, this->key_size_); // DETERMINISTIC FOR TESTING!
}

uchar* CircuitContainer::getKey() {
  return this->key_;
}

uint32 CircuitContainer::getLastCircuitSignNS() {
  return this->last_circuit_sign_time_;
}

uint32 CircuitContainer::getLastCircuitVerifyNS() {
  return this->last_circuit_verify_time_;
}

uint64 CircuitContainer::getLastDirectCallCycles() {
  return this->last_direct_call_cycles_;
}

uint32 CircuitContainer::getLastDirectCallNS() {
  return this->last_direct_call_time_;
}

uint32 CircuitContainer::getViewSizeBits() {
  return (this->num_mul_gates_ * BigIntLib::field_size_bits_);
}

std::string CircuitContainer::getCipherName() {
  return this->cipher_name_;
}

uint32 CircuitContainer::getCipherNumBranches() {
  return this->num_feistel_branches_;
}

uint32 CircuitContainer::getCipherNumMulGates() {
  return this->num_mul_gates_;
}

uint32 CircuitContainer::getCipherNumRounds() {
  return this->num_rounds_;
}

void CircuitContainer::circuitMiMC(word* value_shares_f, word* key_shares_f, uint32 party_size) {
  // Use index 0 ([0]) to indicate first (and only) branch
  word (*shares)[party_size][this->gate_num_words_] = (word (*)[party_size][this->gate_num_words_]) value_shares_f;
  word (*key_shares_field)[party_size][this->gate_num_words_] = (word (*)[party_size][this->gate_num_words_]) key_shares_f;
  // CIRCUIT
  // -------
  // Pre
  word x_squared[party_size][this->gate_num_words_];
  word x_cubed[party_size][this->gate_num_words_];
  memset(x_squared, 0, party_size * this->gate_num_words_ * (WORD_SIZE / 8));
  memset(x_cubed, 0, party_size * this->gate_num_words_ * (WORD_SIZE / 8));
  
  // ADD with key (omit adding constant, since c_0 = 0)
  (this->*add_shared_function_)((word*)shares[0], (word*)key_shares_field[0], (word*)shares[0]);

  // x^3 calculation
  (this->*squ_shared_function_)((word*)shares[0], (word*)shares[0], (word*)x_squared); // x * x = x^2
  (this->*mul_shared_function_)((word*)x_squared, (word*)shares[0], (word*)x_cubed); // x^2 * x = x^3
  this->copyShares((word*)x_cubed, (word*)shares[0], party_size);

  // Rounds
  for(uint32 i = 1; i < this->num_rounds_; i++) {
    // ADD with key
    (this->*add_shared_function_)((word*)shares[0], (word*)key_shares_field, (word*)shares[0]);
  
    // ADD with round constant
    (this->*add_c_shared_function_)((word*)shares[0], this->round_constants_[i], (word*)shares[0]);
  
    // x^3 calculation
    (this->*squ_shared_function_)((word*)shares[0], (word*)shares[0], (word*)x_squared); // x * x = x^2
    (this->*mul_shared_function_)((word*)x_squared, (word*)shares[0], (word*)x_cubed); // x^2 * x = x^3
    this->copyShares((word*)x_cubed, (word*)shares[0], party_size);
  }

  // Post
  // ADD with key
  (this->*add_shared_function_)((word*)shares[0], (word*)key_shares_field, (word*)shares[0]);
  // word temp[this->gate_num_words_];
  // memset(temp, 0, this->gate_num_words_ * (WORD_SIZE / 8));
  // BigIntLib::Add(temp, temp, shares[0]);
  // BigIntLib::Add(temp, temp, shares[1]);
  // BigIntLib::Add(temp, temp, shares[2]);
  // std::cout << "[CIRCUIT] Value: " << BigIntLib::ToString(temp) << std::endl;

  // -------
  for(uint32 i = 0; i < party_size; i++) {
    memcpy(this->value_shares_[i], shares[0][i], this->value_size_);
  }

  // Calculate final shared output value and store it in this->value_ (TEMP, NOT NECESSARY)
  // ---------- DEBUG ----------
  #ifdef VERBOSE
  word tmp_words[this->gate_num_words_];
  std::cout << "[CIRCUIT] Shared output (byte-aligned): ";
  memset(tmp_words, 0, this->gate_size_);
  for(uint32 i = 0; i < party_size; i++) {
    BigIntLib::Add(tmp_words, tmp_words, shares[0][i]);
  }
  std::cout << BigIntLib::ToString(tmp_words);
  std::cout << std::endl;
  #endif
  // ---------- DEBUG ----------
  
  /*
  uint32 offset;
  word temp_words[this->gate_num_words_];
  memset(temp_words, 0, this->value_size_); // gate_size_ == value_size_ for MiMC // TODO!!! IMPORTANT FOR DEBUGGING
  for(uint32 i = 0; i < party_size; i++) {
    BigIntLib::Add(temp_words, temp_words, (word*)(this->value_shares_[i]));
  }
  memcpy(this->value_, temp_words, this->value_size_); // gate_size_ == value_size_ for MiMC
  */
  
}

void CircuitContainer::initMatrixM() {
  // Initialize matrix containing values for the new key schedule

  uint32 matrix_num_words = this->num_feistel_branches_ * this->num_feistel_branches_ * this->gate_num_words_;
  this->matrix_1_ = new word[matrix_num_words];
  memset(this->matrix_1_, 0, matrix_num_words * 8);
  word (*matrix_ptr)[this->num_feistel_branches_][this->gate_num_words_] = (word (*)[this->num_feistel_branches_][this->gate_num_words_]) this->matrix_1_;
  
  for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
    for(uint32 j = 0; j < this->num_feistel_branches_; j++) {
      if(i == j) {
        matrix_ptr[i][j][0] = 2;
      }
      else {
        matrix_ptr[i][j][0] = 1;
      }
    }
  }


  // Print out matrix
  #ifdef VERBOSE
  std::cout << "Matrix M:" << std::endl;
  for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
    for(uint32 j = 0; j < this->num_feistel_branches_; j++) {
      std::cout << BigIntLib::ToString(matrix_ptr[i][j]) << "   ";
    }
    std::cout << std::endl;
  }
  #endif

}

void CircuitContainer::initMatrixMPrime() {
  uint32 matrix_num_words = this->num_feistel_branches_ * this->num_feistel_branches_ * this->gate_num_words_;
  this->matrix_2_ = new word[matrix_num_words];
  memset(this->matrix_2_, 0, matrix_num_words * 8);
  word (*matrix_ptr)[this->num_feistel_branches_][this->gate_num_words_] = (word (*)[this->num_feistel_branches_][this->gate_num_words_]) this->matrix_2_;

  for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
    for(uint32 j = 0; j < this->num_feistel_branches_; j++) {
      if(i == j) {
        matrix_ptr[i][j][0] = 3;
      }
      else {
        matrix_ptr[i][j][0] = 1;
      }
    }
  }

  // Print out matrix
  #ifdef VERBOSE
  std::cout << "Matrix M':" << std::endl;
  for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
    for(uint32 j = 0; j < this->num_feistel_branches_; j++) {
      std::cout << BigIntLib::ToString(matrix_ptr[i][j]) << "   ";
    }
    std::cout << std::endl;
  }
  #endif
}

void CircuitContainer::beforeSign(uchar* x, SignData* sign_data, word* value_shares_f, word* key_shares_f) {

  this->sign_data_ = sign_data;

  // NEW
  // Set this->value_shares_[0] to x, and both this->value_shares_[1] and this->value_shares_[2] to 0
  memcpy(this->value_shares_[0], x, this->value_size_);
  memset(this->value_shares_[1], 0x0, this->value_size_);
  memset(this->value_shares_[2], 0x0, this->value_size_);

  memset(this->key_shares_[0], 0, this->key_size_);
  memset(this->key_shares_[1], 0, this->key_size_);
  memcpy(this->key_shares_[0], (this->sign_data_)->random_tapes_hashs_, this->hash_size_);
  memcpy(this->key_shares_[1], (this->sign_data_)->random_tapes_hashs_ + this->hash_size_, this->hash_size_);

  // Calculate final shares for value and key
  (this->*prepare_shares_field_sign_function_)(x, value_shares_f, key_shares_f);

  // Write x_3 (last key share!) to SignData
  memcpy((this->sign_data_)->x_3_, this->key_shares_[2], this->value_size_);

  #ifdef VERBOSE
  std::cout << "[SIGN] Input private x: ";
  BigIntLib::Print(x, this->value_size_);
  //BigIntLib::Print((word*)this->value_shares_[0], this->value_size_);
  //BigIntLib::Print((word*)this->value_shares_[1], this->value_size_);
  for(uint32 i = 0; i < this->party_size_; i++) {
    std::cout << "[SIGN] Input share " << i << ": ";
    BigIntLib::Print(this->value_shares_[i], this->value_size_);
  }
  #endif

  #ifdef VERBOSE
  for(uint32 i = 0; i < this->party_size_; i++) {
    std::cout << "[SIGN] Key share " << i << ": ";
    BigIntLib::Print(this->key_shares_[i], this->key_size_);
  }
  #endif

  uchar* random_tapes[3];
  random_tapes[0] = (this->sign_data_)->random_tapes_;
  random_tapes[1] = random_tapes[0] + this->random_tape_size_;
  random_tapes[2] = random_tapes[1] + this->random_tape_size_;
  this->prepareRandomNumbers(random_tapes, this->party_size_);

  // Function pointers
  this->add_c_shared_function_ = &CircuitContainer::addCSharedSign;
  this->add_shared_function_ = &CircuitContainer::addSharedSign;
  this->sub_shared_function_ = &CircuitContainer::subSharedSign;
  this->mul_shared_function_ = &CircuitContainer::mulSharedSign;
  this->squ_shared_function_ = &CircuitContainer::squSharedSign;
  // Can be made shorter with one if for first condition
  if(BigIntLib::field_type_ == 1 && BigIntLib::field_size_bits_ == 3) this->squ_shared_experimental_function_ = &CircuitContainer::squSharedExperimental3Sign;
  if(BigIntLib::field_type_ == 1 && BigIntLib::field_size_bits_ == 17) this->squ_shared_experimental_function_ = &CircuitContainer::squSharedExperimental17Sign;
  if(BigIntLib::field_type_ == 1 && BigIntLib::field_size_bits_ == 33) this->squ_shared_experimental_function_ = &CircuitContainer::squSharedExperimental33Sign;
  this->cube_shared_function_ = &CircuitContainer::cubeSharedSign;

  // Reset current mul gate and current intermediate result
  this->current_mul_gate_ = 0;
  this->current_intermediate_result_ = 0;
}

void CircuitContainer::afterSign(word* value_shares_f) {
  // Write shares to bytes
  /*
  this->value_shares_[0] = (this->sign_data_)->y_shares_;
  this->value_shares_[1] = (this->sign_data_)->y_shares_ + this->value_size_;
  this->value_shares_[2] = (this->sign_data_)->y_shares_ + (2 * this->value_size_);
  */
  (this->*output_shares_to_bytes_function_)(value_shares_f, this->party_size_);

  // Write last values to SignData
  for(uint32 i = 0; i < this->party_size_; i++) {
    memcpy((this->sign_data_)->y_shares_ + (i * this->value_size_), this->value_shares_[i], this->value_size_);
    #ifdef VERBOSE
    std::cout << "[SIGN] Output share " << i << ": ";
    BigIntLib::Print((this->sign_data_)->y_shares_ + (i * this->value_size_), this->value_size_);
    #endif
  }
  
  
  // Write y to SignData (only temp, to check values, should all be the same)
  /*
  memcpy((this->sign_data_)->y_, this->value_, this->value_size_);
  #ifdef VERBOSE
  std::cout << "[SIGN] Shared output: ";
  BigIntLib::Print((this->sign_data_)->y_, this->value_size_);
  #endif
  */

  // Destroy unneeded data
  //this->destroyRandomNumbers(random_numbers, this->party_size_);
}

void CircuitContainer::beforeVerify(Proof* p, uchar* x, uchar* y, VerifyData* verify_data, uchar** key_shares, word* value_shares_f, word* key_shares_f, uint32 iteration) {
  this->proof_ = p;
  this->verify_data_ = verify_data;
  this->e_ = p->e_[iteration % this->hash_size_] % this->party_size_; // TODO
  
  this->iteration_ = iteration;

  #ifdef VERBOSE
  std::cout << "[VERIFY] e: " << this->e_ << std::endl;
  #endif

  // Input and the two shares
  // Remark: shares[0] and shares[1] are BOTH needed and should be written throughout the computation! each given View value after a MUL must be written to shares[1]!

  memset(this->key_shares_[0], 0, this->key_size_);
  memset(this->key_shares_[1], 0, this->key_size_);
  memset(this->value_shares_[0], 0, this->value_size_);
  memset(this->value_shares_[1], 0, this->value_size_);

  if(this->e_ == 0) {
    memcpy(this->key_shares_[0], p->zs_[iteration]->k_1_hash_, this->hash_size_);
    memcpy(this->key_shares_[1], p->zs_[iteration]->k_2_hash_, this->hash_size_);
    memcpy(this->value_shares_[0], x, this->value_size_);
    /*
    // Generic solution for security < 256 bits (a bit slower due to std::min, could be stored in var at the beginning)
    memcpy(this->value_shares_[0], p->zs_[iteration]->k_1_hash_, std::min(this->value_size_, this->hash_size_));
    memcpy(this->value_shares_[1], p->zs_[iteration]->k_2_hash_, std::min(this->value_size_, this->hash_size_));
    */
    
  }
  else if(this->e_ == 1) {
    memcpy(this->key_shares_[0], p->zs_[iteration]->k_2_hash_, this->hash_size_);
    memcpy(this->key_shares_[1], p->zs_[iteration]->x_3_, this->value_size_);
    /*
    // Generic solution for security < 256 bits (a bit slower due to std::min, could be stored in var at the beginning)
    memcpy(this->value_shares_[0], p->zs_[iteration]->k_2_hash_, std::min(this->value_size_, this->hash_size_));
    memcpy(this->value_shares_[1], p->zs_[iteration]->x_3_, this->value_size_);
    */
    
  }
  else if(this->e_ == 2) {
    memcpy(this->key_shares_[0], p->zs_[iteration]->x_3_, this->value_size_);
    memcpy(this->key_shares_[1], p->zs_[iteration]->k_1_hash_, this->hash_size_);
    memcpy(this->value_shares_[1], x, this->value_size_);
    /*
    // Generic solution for security < 256 bits (a bit slower due to std::min, could be stored in var at the beginning)
    memcpy(this->value_shares_[0], p->zs_[iteration]->x_3_, this->value_size_);
    memcpy(this->value_shares_[1], p->zs_[iteration]->k_1_hash_, std::min(this->value_size_, this->hash_size_));
    */
  }
  else {
    std::cout << "Error: e not in {0, 1, 2}" << std::endl;
  }

  #ifdef VERBOSE
  for(uint32 i = 0; i < this->party_size_ - 1; i++) {
    std::cout << "[VERIFY] Input share " << i << ": ";
    BigIntLib::Print(this->value_shares_[i], this->value_size_);
  }
  #endif

  #ifdef VERBOSE
  for(uint32 i = 0; i < this->party_size_ - 1; i++) {
    std::cout << "[VERIFY] Key share " << i << ": ";
    BigIntLib::Print(this->key_shares_[i], this->key_size_);
  }
  #endif

  uchar* random_tapes[2];
  random_tapes[0] = p->zs_[iteration]->k_1_;
  random_tapes[1] = p->zs_[iteration]->k_2_;
  this->prepareRandomNumbers(random_tapes, this->party_size_ - 1);

  (this->*prepare_shares_field_verify_function_)(value_shares_f, key_shares_f);

  // Function pointers
  this->add_c_shared_function_ = &CircuitContainer::addCSharedVerify;
  this->add_shared_function_ = &CircuitContainer::addSharedVerify;
  this->sub_shared_function_ = &CircuitContainer::subSharedVerify;
  this->mul_shared_function_ = &CircuitContainer::mulSharedVerify;
  this->squ_shared_function_ = &CircuitContainer::squSharedVerify;
  // Can be made shorter with one if for first condition
  if(BigIntLib::field_type_ == 1 && BigIntLib::field_size_bits_ == 3) this->squ_shared_experimental_function_ = &CircuitContainer::squSharedExperimental3Verify;
  if(BigIntLib::field_type_ == 1 && BigIntLib::field_size_bits_ == 17) this->squ_shared_experimental_function_ = &CircuitContainer::squSharedExperimental17Verify;
  if(BigIntLib::field_type_ == 1 && BigIntLib::field_size_bits_ == 33) this->squ_shared_experimental_function_ = &CircuitContainer::squSharedExperimental33Verify;
  this->cube_shared_function_ = &CircuitContainer::cubeSharedVerify;

  // Reset current mul gate and current intermediate result
  this->current_mul_gate_ = 0;
  this->current_intermediate_result_ = 0;
}

void CircuitContainer::afterVerify(uchar* y, word* value_shares_f) {
  // Write shares to bytes
  (this->*output_shares_to_bytes_function_)(value_shares_f, this->party_size_ - 1);

  memcpy((this->verify_data_)->y_share_, this->value_shares_[0], this->value_size_);

  #ifdef VERBOSE
  std::cout << "[VERIFY] Output: ";
  BigIntLib::Print(y, this->value_size_);
  std::cout << "[VERIFY] Calculated Output share: ";
  BigIntLib::Print((this->verify_data_)->y_share_, this->value_size_);
  std::cout << "[VERIFY] Given Output share: ";
  BigIntLib::Print(((this->proof_)->zs_[this->iteration_])->y_share_, this->value_size_);
  #endif

  // Calculate last share
  (this->*verify_calc_last_share_function_)(y, value_shares_f);

  #ifdef VERBOSE
  std::cout << "[VERIFY] Calculated Final Output share: ";
  BigIntLib::Print((this->verify_data_)->y_e2_, this->value_size_);
  #endif

  // Destroy unneeded data
  //this->destroyRandomNumbers(random_numbers, this->party_size_ - 1);
}

void CircuitContainer::prepareRoundConstants(uint32 length) {
  if(this->round_constants_ != NULL) {
    for(uint32 i = 0; i < length; i++) {
      delete[] this->round_constants_[i];
    }
    delete[] this->round_constants_;
  }
  this->round_constants_ = new word*[length];
  uint32 msb_word_index = this->gate_num_words_ - 1;
  for(uint32 i = 0; i < length; i++) {
    this->round_constants_[i] = new word[this->gate_num_words_];
    //BigIntLib::FillRandom(this->round_constants_[i], this->gate_size_);
    //memset(this->round_constants_[i], 0, (this->gate_num_words_) * (WORD_SIZE / 8));
    RAND_bytes((uchar*)this->round_constants_[i], this->gate_size_);
    //memset(this->round_constants_[i], 0xCD, this->gate_size_); // DETERMINISTIC FOR TESING!
    // Keep in field (this takes a lot of time...)
    this->round_constants_[i][msb_word_index] &= BigIntLib::msb_word_mask_;
    BigIntLib::TryReduce(this->round_constants_[i]); // Round constants are public and fixed
  }
}

void CircuitContainer::prepareRandomNumbers(uchar** random_tapes, uint32 party_size) {
  //uint32 offset;
  //word* tmp_words;
  //uint32 msb_word_index = this->gate_num_words_ - 1;
  for(uint32 i = 0; i < party_size; i++) {
    // Use random tape of this party as seed and set random numbers for mul gates
    BigIntLib::GetRandomFieldElements(this->random_numbers_[i], random_tapes[i], this->num_mul_gates_);
  }
}

void CircuitContainer::destroyRandomNumbers() {
  for(uint32 i = 0; i < this->party_size_; i++) {
    delete[] this->random_numbers_[i];
  }
  delete[] this->random_numbers_;
}

void CircuitContainer::generateNewSubkeys(word* key_shares_f, word* key_share_f_temp, uint32 num_rounds_remaining) {
  word (*matrix_ptr)[this->num_feistel_branches_][this->gate_num_words_] = (word (*)[this->num_feistel_branches_][this->gate_num_words_]) this->matrix_1_;
  word (*key_shares_field)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) key_shares_f;
  word (*key_shares_field_temp)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) key_share_f_temp;
  memset(key_shares_field_temp, 0, this->num_feistel_branches_ * this->gate_size_);

  // Cache-friendly matrix-vector multiplication, can probably be made faster by splitting into sections
  // Cache line size is typically 64 bytes, so maybe do it column-wise and with floor(512 / (this->gate_size_ * 8)) columns of M in each step
  // Remark: M is in the heap!
  /*
  word temp[this->gate_num_words_];
  uint32 num_subkeys_needed = std::min(this->num_feistel_branches_, num_rounds_remaining);
  for(uint32 i = 0; i < num_subkeys_needed; i++) {
    BigIntLib::Mul(key_shares_field_temp[i], matrix_ptr[i][0], key_shares_field[i]);
  }
  for(uint32 j = 1; j < this->num_feistel_branches_; j++) {
    for(uint32 i = 0; i < num_subkeys_needed; i++) {
      BigIntLib::Mul(temp, matrix_ptr[i][j], key_shares_field[i]);
      BigIntLib::Add(key_shares_field_temp[i], key_shares_field_temp[i], temp);
    }
  }
  */

  // Strictly diagonally-dominant matrix with M_{i,j} = 16 for i = j and M_{i,j} = 1 for i != j (16-bit prime field test)
  word temp[this->gate_num_words_];
  memset(temp, 0x0, this->gate_size_);
  uint32 num_subkeys_needed = std::min(this->num_feistel_branches_, num_rounds_remaining);
  for(uint32 j = 0; j < this->num_feistel_branches_; j++) {
    for(uint32 i = 0; i < num_subkeys_needed; i++) {
      /*
      if(j == i) {
        BigIntLib::Add(key_shares_field_temp[i], key_shares_field_temp[i], key_shares_field[i]);
        BigIntLib::AddSpec(key_shares_field_temp[i], key_shares_field_temp[i], key_shares_field[i]);
        //BigIntLib::MulSpec(temp, key_shares_field[i]);
        //BigIntLib::AddSpec(key_shares_field_temp[i], key_shares_field_temp[i], temp);
      }
      else {
        BigIntLib::AddSpec(key_shares_field_temp[i], key_shares_field_temp[i], key_shares_field[i]);
      }
      */
      if(i == j) { // if(matrix_ptr[i][j][0] == 2) {
        //BigIntLib::Add(key_shares_field_temp[i], key_shares_field_temp[i], key_shares_field[j]);
        memset(temp, 0x0, this->gate_size_);
        BigIntLib::Times2(temp, key_shares_field[j]);
        BigIntLib::Add(key_shares_field_temp[i], key_shares_field_temp[i], temp);
      }
      else {
        BigIntLib::Add(key_shares_field_temp[i], key_shares_field_temp[i], key_shares_field[j]);
      }
    }
  }
  // Reduce
  /*
  for(uint32 i = 0; i < num_subkeys_needed; i++) {
    BigIntLib::SolinasReduc(key_shares_field_temp[i], key_shares_field_temp[i]);
  }
  */

  // Copy and reset
  memcpy(key_shares_field, key_shares_field_temp, this->num_feistel_branches_ * this->gate_size_);
  memset(key_shares_field_temp, 0x0, this->num_feistel_branches_ * this->gate_size_);
}

void CircuitContainer::generateNewSubkeysSharedPF(word* key_shares_f, word* key_share_f_temp, uint32 num_rounds_remaining, uint32 party_size) {
  word (*matrix_ptr)[this->num_feistel_branches_][this->gate_num_words_] = (word (*)[this->num_feistel_branches_][this->gate_num_words_]) this->matrix_1_;
  word (*key_shares_field)[party_size][this->gate_num_words_] = (word (*)[party_size][this->gate_num_words_]) key_shares_f;
  word (*key_shares_field_temp)[party_size][this->gate_num_words_] = (word (*)[party_size][this->gate_num_words_]) key_share_f_temp;
  memset(key_shares_field_temp, 0, this->num_feistel_branches_ * party_size * this->gate_size_);
  
  word temp[this->gate_num_words_];
  word temp_2[this->gate_num_words_];
  memset(temp, 0x0, this->gate_size_);
  memset(temp_2, 0x0, this->gate_size_);
  for(uint32 k = 0; k < party_size; k++) {
    memset(temp, 0x0, this->gate_size_);
    for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
      BigIntLib::Add(temp, temp, key_shares_field[i][k]);
    }
    // Add/Sub S_x to/from each state, do correcting computations in binary field
    //memset(temp2, 0x0, this->gate_size_);
    for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
      BigIntLib::Add(key_shares_field[i][k], key_shares_field[i][k], temp);
    }
  }
}

void CircuitContainer::generateNewSubkeysSharedBF(word* key_shares_f, word* key_share_f_temp, uint32 num_rounds_remaining, uint32 party_size) {
  word (*matrix_ptr)[this->num_feistel_branches_][this->gate_num_words_] = (word (*)[this->num_feistel_branches_][this->gate_num_words_]) this->matrix_1_;
  word (*key_shares_field)[party_size][this->gate_num_words_] = (word (*)[party_size][this->gate_num_words_]) key_shares_f;
  word (*key_shares_field_temp)[party_size][this->gate_num_words_] = (word (*)[party_size][this->gate_num_words_]) key_share_f_temp;
  memset(key_shares_field_temp, 0, this->num_feistel_branches_ * party_size * this->gate_size_);
  
  word temp[this->gate_num_words_];
  word temp_2[this->gate_num_words_];
  memset(temp, 0x0, this->gate_size_);
  memset(temp_2, 0x0, this->gate_size_);
  for(uint32 k = 0; k < party_size; k++) {
    memset(temp, 0x0, this->gate_size_);
    for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
      BigIntLib::Add(temp, temp, key_shares_field[i][k]);
    }
    // Add/Sub S_x to/from each state, do correcting computations in binary field
    //memset(temp2, 0x0, this->gate_size_);
    for(uint32 i = 0; i < this->num_feistel_branches_; i++) { 
      BigIntLib::Times2(temp_2, key_shares_field[i][k]);
      BigIntLib::Sub(key_shares_field[i][k], key_shares_field[i][k], temp);
      BigIntLib::Add(key_shares_field[i][k], key_shares_field[i][k], temp_2);
    }
  }
}

// word* shares contains the words of ALL shares (share_party_1_word_1, share_party_1_word_2, ... ,share_party_n_word_m-1, share_party_n_word_m) for n parties and m words
void CircuitContainer::addSharedSign(word* a_shares, word* b_shares, word* c_shares) {
  
  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*b_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) b_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;
  BigIntLib::Add(c_shares_p[0], a_shares_p[0], b_shares_p[0]);
  BigIntLib::Add(c_shares_p[1], a_shares_p[1], b_shares_p[1]);
  BigIntLib::Add(c_shares_p[2], a_shares_p[2], b_shares_p[2]);
}

// word* shares contains the words of ALL shares (share_party_1_word_1, share_party_1_word_2, ... ,share_party_n_word_m-1, share_party_n_word_m) for n parties and m words
void CircuitContainer::addSharedVerify(word* a_shares, word* b_shares, word* c_shares) {
  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*b_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) b_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;
  BigIntLib::Add(c_shares_p[0], a_shares_p[0], b_shares_p[0]);
  BigIntLib::Add(c_shares_p[1], a_shares_p[1], b_shares_p[1]);
}

void CircuitContainer::subSharedSign(word* a_shares, word* b_shares, word* c_shares) {
  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*b_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) b_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;
  BigIntLib::Sub(c_shares_p[0], a_shares_p[0], b_shares_p[0]);
  BigIntLib::Sub(c_shares_p[1], a_shares_p[1], b_shares_p[1]);
  BigIntLib::Sub(c_shares_p[2], a_shares_p[2], b_shares_p[2]);
}

void CircuitContainer::subSharedVerify(word* a_shares, word* b_shares, word* c_shares) {
  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*b_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) b_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;
  BigIntLib::Sub(c_shares_p[0], a_shares_p[0], b_shares_p[0]);
  BigIntLib::Sub(c_shares_p[1], a_shares_p[1], b_shares_p[1]);
}

// word* shares contains the words of ALL shares (share_party_1_word_1, share_party_1_word_2, ... ,share_party_n_word_m-1, share_party_n_word_m) for n parties and m words
void CircuitContainer::addCSharedSign(word* a_shares, word* b, word* c_shares) {
  uint32 offset;
  BigIntLib::Add(c_shares, a_shares, b);
}

// word* shares contains the words of ALL shares (share_party_1_word_1, share_party_1_word_2, ... ,share_party_n_word_m-1, share_party_n_word_m) for n parties and m words
void CircuitContainer::addCSharedVerify(word* a_shares, word* b, word* c_shares) {
  if(this->e_ == 0) {
    BigIntLib::Add(c_shares, a_shares, b);
    //memcpy(c_shares + this->gate_num_words_, a_shares + this->gate_num_words_, this->gate_size_);
  }
  //else if(this->e_ == 1) {
    //memcpy(c_shares, a_shares, this->gate_size_);
    //memcpy(c_shares + this->gate_num_words_, a_shares + this->gate_num_words_, this->gate_size_);
  //}
  else if(this->e_ == 2) {
    //memcpy(c_shares, a_shares, this->gate_size_);
    BigIntLib::Add(c_shares + this->gate_num_words_, a_shares + this->gate_num_words_, b);
  }
}

void CircuitContainer::squSharedSign(word* a_shares, word* b_shares, word* c_shares) {

  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*b_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) b_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;

  word temp[this->gate_num_words_];
  uint32 mul_gate_random_pointer = this->current_mul_gate_ * this->gate_num_words_;
  uint32 mul_gate_pointer_uchar = this->current_mul_gate_ * this->gate_size_;
  
  // Share 1
  BigIntLib::Mul(temp, a_shares_p[0], b_shares_p[0]);
  BigIntLib::Mul(c_shares_p[0], a_shares_p[1], b_shares_p[0]);
  BigIntLib::Add(temp, temp, c_shares_p[0]);
  BigIntLib::Add(temp, temp, c_shares_p[0]);
  BigIntLib::Add(temp, temp, this->random_numbers_[0] + mul_gate_random_pointer);
  BigIntLib::Sub(c_shares_p[0], temp, this->random_numbers_[1] + mul_gate_random_pointer);
  // Write c share to SignData
  memcpy((this->sign_data_)->views_[0] + mul_gate_pointer_uchar, c_shares_p[0], BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  //memcpy(view_pointer[0][mul_gate_random_pointer], c_shares_p[0], BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_

  // Share 2
  BigIntLib::Mul(temp, a_shares_p[1], b_shares_p[1]);
  BigIntLib::Mul(c_shares_p[1], a_shares_p[2], b_shares_p[1]);
  BigIntLib::Add(temp, temp, c_shares_p[1]);
  BigIntLib::Add(temp, temp, c_shares_p[1]);
  BigIntLib::Add(temp, temp, this->random_numbers_[1] + mul_gate_random_pointer);
  BigIntLib::Sub(c_shares_p[1], temp, this->random_numbers_[2] + mul_gate_random_pointer);
  // Write c share to SignData
  memcpy((this->sign_data_)->views_[1] + mul_gate_pointer_uchar, c_shares_p[1], BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  //memcpy(view_pointer[1][mul_gate_random_pointer], c_shares_p[1], BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_

  // --- New ZKB++ optimization
  BigIntLib::Sub(temp, this->intermediate_results_[this->current_intermediate_result_], c_shares_p[0]);
  BigIntLib::Sub(c_shares_p[2], temp, c_shares_p[1]);
  memcpy((this->sign_data_)->views_[2] + mul_gate_pointer_uchar, c_shares_p[2], BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  this->current_intermediate_result_++;
  // ---

  this->current_mul_gate_++;
}

void CircuitContainer::squSharedExperimental3Sign(word* a_shares, word* b_shares, word* c_shares) {

  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;

  // Procedure for each party separately:
  // 1. Apply linear squaring algorithm
  // 2. Apply steps for modular reduction

  doubleword r;
  word u_0, u_1, u_2, T;
  word c0, c1;
  // Party 1
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[0][0] & 0xFF;
  r = this->squaring_precomp_[u_0];
  // Reduction
  // p(x) = x^3 + x + 1
  c0 = r & 0x7; // LS 3 bits
  c1 = (r & 0x18) >> 3; // MS 2 bits
  // Add c1 to bits 1 and 0 of c0 (c1 doesn't neet to be done first, because it's not affected)
  c_shares_p[0][0] = c0 ^ (c1 << 1) ^ c1;
  
  // Party 2
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[1][0] & 0xFF;
  r = this->squaring_precomp_[u_0];
  // Reduction
  // p(x) = x^3 + x + 1
  c0 = r & 0x7; // LS 3 bits
  c1 = (r & 0x18) >> 3; // MS 2 bits
  // Add c1 to bits 1 and 0 of c0 (c1 doesn't neet to be done first, because it's not affected)
  c_shares_p[1][0] = c0 ^ (c1 << 1) ^ c1;

  
  // Party 3
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[2][0] & 0xFF;
  r = this->squaring_precomp_[u_0];
  // Reduction
  // p(x) = x^3 + x + 1
  c0 = r & 0x7; // LS 3 bits
  c1 = (r & 0x18) >> 3; // MS 2 bits
  // Add c1 to bits 1 and 0 of c0 (c1 doesn't neet to be done first, because it's not affected)
  c_shares_p[2][0] = c0 ^ (c1 << 1) ^ c1;
  

  // --- Necessary for new ZKB++ optimization (here, the old method tends to be faster for squaring)
  this->current_intermediate_result_++;
  // ---
}

void CircuitContainer::squSharedExperimental17Sign(word* a_shares, word* b_shares, word* c_shares) {

  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;

  // Procedure for each party separately:
  // 1. Apply linear squaring algorithm
  // 2. Apply steps for modular reduction

  doubleword r, t;
  word u_0, u_1, u_2, T;
  word r_0_n_1, r_n_deg_r;
  word x_m_1 = 0x9; // = 2^3 + 1
  word c0, c1;
  // Party 1
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[0][0] & 0xFF;
  u_1 = (a_shares_p[0][0] & 0xFF00) >> 8;
  u_2 = (a_shares_p[0][0] & 0xFF0000) >> 16;
  r = this->squaring_precomp_[u_2];
  r = (r << 32) | (this->squaring_precomp_[u_1] << 16) | this->squaring_precomp_[u_0];
  word r_tmp = r;
  // Reduction
  // p(x) = x^17 + x^3 + 1
  r_0_n_1 = r & 0x1FFFF;
  r_n_deg_r = r >> 17;
  t = 0;
  asm("pclmulqdq %2, %1, %0;"
    : "=x"(t)
    : "x"(r_n_deg_r), "i"(0), "0"(x_m_1)
    );
  r = r_0_n_1 ^ t;
  t = r >> 17;
  r = r ^ (t << 3) ^ t;
  c_shares_p[0][0] = r & 0x1FFFF;
  
  /*
  c0 = r_tmp & 0x1FFFF; // LS 17 bits
  c1 = r_tmp >> 17; // MS 16 bits
  
  c1 = c1 ^ (c1 >> 14);
  // Add c1 to bits 3 and 0 of c0
  c0 = c0 ^ (c1 << 3) ^ c1;
  
  // Build result
  c_shares_p[0][0] = c0 & 0x1FFFF;
  */

  // Party 2
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[1][0] & 0xFF;
  u_1 = (a_shares_p[1][0] & 0xFF00) >> 8;
  u_2 = (a_shares_p[1][0] & 0xFF0000) >> 16;
  r = this->squaring_precomp_[u_2];
  r = (r << 32) | (this->squaring_precomp_[u_1] << 16) | this->squaring_precomp_[u_0];
  // Reduction
  // p(x) = x^17 + x^3 + 1
  
  r_0_n_1 = r & 0x1FFFF;
  r_n_deg_r = r >> 17;
  t = 0;
  asm("pclmulqdq %2, %1, %0;"
    : "=x"(t)
    : "x"(r_n_deg_r), "i"(0), "0"(x_m_1)
    );
  r = r_0_n_1 ^ t;
  t = r >> 17;
  r = r ^ (t << 3) ^ t;
  c_shares_p[1][0] = r & 0x1FFFF;
  
  /*
  c0 = r & 0x1FFFF; // LS 17 bits
  c1 = r >> 17; // MS 16 bits
  
  c1 = c1 ^ (c1 >> 14);
  // Add c1 to bits 3 and 0 of c0
  c0 = c0 ^ (c1 << 3) ^ c1;
  
  // Build result
  c_shares_p[1][0] = c0 & 0x1FFFF;
  */

  /*
  // Party 3
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[2][0] & 0xFF;
  u_1 = (a_shares_p[2][0] & 0xFF00) >> 8;
  u_2 = (a_shares_p[2][0] & 0xFF0000) >> 16;
  r = this->squaring_precomp_[u_2];
  r = (r << 32) | (this->squaring_precomp_[u_1] << 16) | this->squaring_precomp_[u_0];
  // Reduction
  // p(x) = x^17 + x^3 + 1
  
  r_0_n_1 = r & 0x1FFFF;
  r_n_deg_r = r >> 17;
  t = 0;
  asm("pclmulqdq %2, %1, %0;"
    : "=x"(t)
    : "x"(r_n_deg_r), "i"(0), "0"(x_m_1)
    );
  r = r_0_n_1 ^ t;
  t = r >> 17;
  r = r ^ (t << 3) ^ t;
  c_shares_p[2][0] = r & 0x1FFFF;
  */

  // --- New ZKB++ optimization
  word temp[this->gate_num_words_];
  BigIntLib::Sub(temp, this->intermediate_results_[this->current_intermediate_result_], c_shares_p[0]);
  BigIntLib::Sub(c_shares_p[2], temp, c_shares_p[1]);
  this->current_intermediate_result_++;
  // ---
  
  /*
  c0 = r & 0x1FFFF; // LS 17 bits
  c1 = r >> 17; // MS 16 bits
  
  c1 = c1 ^ (c1 >> 14);
  // Add c1 to bits 3 and 0 of c0
  c0 = c0 ^ (c1 << 3) ^ c1;
  
  // Build result
  c_shares_p[2][0] = c0 & 0x1FFFF;
  */
}

void CircuitContainer::squSharedExperimental33Sign(word* a_shares, word* b_shares, word* c_shares) {

  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;

  // Procedure for each party separately:
  // 1. Apply linear squaring algorithm
  // 2. Apply steps for modular reduction

  doubleword r;
  word u_0, u_1, u_2, u_3, c0, c1, T;
  // Party 1
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[0][0] & 0xFF;
  u_1 = (a_shares_p[0][0] & 0xFF00) >> 8;
  u_2 = (a_shares_p[0][0] & 0xFF0000) >> 16;
  u_3 = (a_shares_p[0][0] & 0xFF000000) >> 24;
  r = (this->squaring_precomp_[u_3] << 16) | this->squaring_precomp_[u_2];
  r = (r << 32) | (this->squaring_precomp_[u_1] << 16) | this->squaring_precomp_[u_0];
  // Last 1 bit (MS 33th bit)
  u_0 = (a_shares_p[0][0] & 0x100000000) >> 32;
  r = r | ((doubleword)(this->squaring_precomp_[u_0]) << 64);
  // Reduction
  c0 = r & 0x1FFFFFFFF; // LS 33 bits
  c1 = r >> 33; // MS 32 bits
  T = c1;
  c1 = c1 ^ (T >> 27) ^ (T >> 30); // 27 = 33 - 6, 30 = 33 - 3, 32 = 33 - 1 (omitted, all zeros), x^0 does not affect c1
  T = c1;
  c_shares_p[0][0] = c0 ^ ((T << 6) & 0x1FFFFFFFF) ^ ((T << 3) & 0x1FFFFFFFF) ^ ((T << 1) & 0x1FFFFFFFF) ^ T; // = c0, for x^6, x^3, x^1, x^0
  

  // Party 2
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[1][0] & 0xFF;
  u_1 = (a_shares_p[1][0] & 0xFF00) >> 8;
  u_2 = (a_shares_p[1][0] & 0xFF0000) >> 16;
  u_3 = (a_shares_p[1][0] & 0xFF000000) >> 24;
  r = (this->squaring_precomp_[u_3] << 16) | this->squaring_precomp_[u_2];
  r = (r << 32) | (this->squaring_precomp_[u_1] << 16) | this->squaring_precomp_[u_0];
  // Last 1 bit (MS 33th bit)
  u_0 = (a_shares_p[1][0] & 0x100000000) >> 32;
  r = r | ((doubleword)(this->squaring_precomp_[u_0]) << 64);
  // Reduction
  c0 = r & 0x1FFFFFFFF; // LS 33 bits
  c1 = r >> 33; // MS 32 bits
  T = c1;
  c1 = c1 ^ (T >> 27) ^ (T >> 30); // 27 = 33 - 6, 30 = 33 - 3, 32 = 33 - 1 (omitted, all zeros), x^0 does not affect c1
  T = c1;
  c_shares_p[1][0] = c0 ^ ((T << 6) & 0x1FFFFFFFF) ^ ((T << 3) & 0x1FFFFFFFF) ^ ((T << 1) & 0x1FFFFFFFF) ^ T; // = c0, for x^6, x^3, x^1, x^0

  /*
  // Party 3
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[2][0] & 0xFF;
  u_1 = (a_shares_p[2][0] & 0xFF00) >> 8;
  u_2 = (a_shares_p[2][0] & 0xFF0000) >> 16;
  u_3 = (a_shares_p[2][0] & 0xFF000000) >> 24;
  r = (this->squaring_precomp_[u_3] << 16) | this->squaring_precomp_[u_2];
  r = (r << 32) | (this->squaring_precomp_[u_1] << 16) | this->squaring_precomp_[u_0];
  // Last 1 bit (MS 33th bit)
  u_0 = (a_shares_p[2][0] & 0x100000000) >> 32;
  r = r | ((doubleword)(this->squaring_precomp_[u_0]) << 64);
  // Reduction
  c0 = r & 0x1FFFFFFFF; // LS 33 bits
  c1 = r >> 33; // MS 32 bits
  T = c1;
  c1 = c1 ^ (T >> 27) ^ (T >> 30); // 27 = 33 - 6, 30 = 33 - 3, 32 = 33 - 1 (omitted, all zeros), x^0 does not affect c1
  T = c1;
  c_shares_p[2][0] = c0 ^ ((T << 6) & 0x1FFFFFFFF) ^ ((T << 3) & 0x1FFFFFFFF) ^ ((T << 1) & 0x1FFFFFFFF) ^ T; // = c0, for x^6, x^3, x^1, x^0
  */

  // --- New ZKB++ optimization
  word temp[this->gate_num_words_];
  BigIntLib::Sub(temp, this->intermediate_results_[this->current_intermediate_result_], c_shares_p[0]);
  BigIntLib::Sub(c_shares_p[2], temp, c_shares_p[1]);
  this->current_intermediate_result_++;
  // ---

}

// word* shares contains the words of ALL shares (share_party_1_word_1, share_party_1_word_2, ... ,share_party_n_word_m-1, share_party_n_word_m) for n parties and m words
void CircuitContainer::mulSharedSign(word* a_shares, word* b_shares, word* c_shares) {

  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*b_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) b_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;

  word temp[this->gate_num_words_];
  // uint32 i_next;
  // uint32 offset_1;
  // uint32 offset_2;
  uint32 mul_gate_random_pointer = this->current_mul_gate_ * this->gate_num_words_;
  uint32 mul_gate_pointer_uchar = this->current_mul_gate_ * this->gate_size_;
  
  // Share 1
  BigIntLib::Mul(temp, a_shares_p[0], b_shares_p[0]);
  BigIntLib::Mul(c_shares_p[0], a_shares_p[1], b_shares_p[0]);
  BigIntLib::Add(temp, temp, c_shares_p[0]);
  BigIntLib::Mul(c_shares_p[0], a_shares_p[0], b_shares_p[1]);
  BigIntLib::Add(temp, temp, c_shares_p[0]);
  BigIntLib::Add(temp, temp, this->random_numbers_[0] + mul_gate_random_pointer);
  BigIntLib::Sub(c_shares_p[0], temp, this->random_numbers_[1] + mul_gate_random_pointer);
  // Write c share to SignData
  memcpy((this->sign_data_)->views_[0] + mul_gate_pointer_uchar, c_shares_p[0], BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  //memcpy(view_pointer[0][mul_gate_random_pointer], c_shares_p[0], BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_

  // Share 2
  BigIntLib::Mul(temp, a_shares_p[1], b_shares_p[1]);
  BigIntLib::Mul(c_shares_p[1], a_shares_p[2], b_shares_p[1]);
  BigIntLib::Add(temp, temp, c_shares_p[1]);
  BigIntLib::Mul(c_shares_p[1], a_shares_p[1], b_shares_p[2]);
  BigIntLib::Add(temp, temp, c_shares_p[1]);
  BigIntLib::Add(temp, temp, this->random_numbers_[1] + mul_gate_random_pointer);
  BigIntLib::Sub(c_shares_p[1], temp, this->random_numbers_[2] + mul_gate_random_pointer);
  // Write c share to SignData
  memcpy((this->sign_data_)->views_[1] + mul_gate_pointer_uchar, c_shares_p[1], BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  //memcpy(view_pointer[1][mul_gate_random_pointer], c_shares_p[1], BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_

  // --- New ZKB++ optimization
  BigIntLib::Sub(temp, this->intermediate_results_[this->current_intermediate_result_], c_shares_p[0]);
  BigIntLib::Sub(c_shares_p[2], temp, c_shares_p[1]);
  memcpy((this->sign_data_)->views_[2] + mul_gate_pointer_uchar, c_shares_p[2], BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  this->current_intermediate_result_++;
  // ---

  this->current_mul_gate_++;
}

void CircuitContainer::squSharedVerify(word* a_shares, word* b_shares, word* c_shares) {
  word temp[this->gate_num_words_];
  uint32 offset = this->gate_num_words_;
  uint32 mul_gate_random_pointer = this->current_mul_gate_ * this->gate_num_words_;
  uint32 mul_gate_pointer_uchar = this->current_mul_gate_ * this->gate_size_;
  BigIntLib::Mul(temp, a_shares, b_shares);
  BigIntLib::Mul(c_shares, a_shares + offset, b_shares);
  BigIntLib::Add(temp, temp, c_shares);
  BigIntLib::Add(temp, temp, c_shares);
  BigIntLib::Add(temp, temp, this->random_numbers_[0] + mul_gate_random_pointer);
  BigIntLib::Sub(c_shares, temp, this->random_numbers_[1] + mul_gate_random_pointer);
  // Get value from view and store in c_shares[1]
  memcpy(c_shares + offset, this->proof_->zs_[this->iteration_]->view_ + mul_gate_pointer_uchar, BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  // Write c_shares[0] to VerifyData (this is part of the View computed now)
  memcpy((this->verify_data_)->view_ + mul_gate_pointer_uchar, c_shares, BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  this->current_mul_gate_++;
}

void CircuitContainer::squSharedExperimental3Verify(word* a_shares, word* b_shares, word* c_shares) {

  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;

  // Procedure for each party separately:
  // 1. Apply linear squaring algorithm
  // 2. Apply steps for modular reduction

  doubleword r;
  word u_0, u_1, u_2, T;
  word c0, c1;
  // Party 1
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[0][0] & 0xFF;
  r = this->squaring_precomp_[u_0];
  // Reduction
  // p(x) = x^3 + x + 1
  c0 = r & 0x7; // LS 3 bits
  c1 = (r & 0x18) >> 3; // MS 2 bits
  // Add c1 to bits 1 and 0 of c0 (c1 doesn't neet to be done first, because it's not affected)
  c_shares_p[0][0] = c0 ^ (c1 << 1) ^ c1;
  
  // Party 2
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[1][0] & 0xFF;
  r = this->squaring_precomp_[u_0];
  // Reduction
  // p(x) = x^3 + x + 1
  c0 = r & 0x7; // LS 3 bits
  c1 = (r & 0x18) >> 3; // MS 2 bits
  // Add c1 to bits 1 and 0 of c0 (c1 doesn't neet to be done first, because it's not affected)
  c_shares_p[1][0] = c0 ^ (c1 << 1) ^ c1;

}

void CircuitContainer::squSharedExperimental17Verify(word* a_shares, word* b_shares, word* c_shares) {

  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;

  // Procedure for each party separately:
  // 1. Apply linear squaring algorithm
  // 2. Apply steps for modular reduction

  doubleword r;
  word u_0, u_1, u_2, T;
  // Party 1
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[0][0] & 0xFF;
  u_1 = (a_shares_p[0][0] & 0xFF00) >> 8;
  u_2 = (a_shares_p[0][0] & 0xFF0000) >> 16;
  r = this->squaring_precomp_[u_2];
  r = (r << 32) | (this->squaring_precomp_[u_1] << 16) | this->squaring_precomp_[u_0];
  // Reduction
  // p(x) = x^17 + x^3 + 1
  word r_0_n_1 = r & 0x1FFFF;
  word r_n_deg_r = r >> 17;
  word x_m_1 = 0x9; // = 2^3 + 1
  doubleword t = 0;
  asm("pclmulqdq %2, %1, %0;"
    : "=x"(t)
    : "x"(r_n_deg_r), "i"(0), "0"(x_m_1)
    );
  r = r_0_n_1 ^ t;
  t = r >> 17;
  r = r ^ (t << 3) ^ t;
  c_shares_p[0][0] = r & 0x1FFFF;

  // Party 2
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[1][0] & 0xFF;
  u_1 = (a_shares_p[1][0] & 0xFF00) >> 8;
  u_2 = (a_shares_p[1][0] & 0xFF0000) >> 16;
  r = this->squaring_precomp_[u_2];
  r = (r << 32) | (this->squaring_precomp_[u_1] << 16) | this->squaring_precomp_[u_0];
  // Reduction
  // p(x) = x^17 + x^3 + 1
  r_0_n_1 = r & 0x1FFFF;
  r_n_deg_r = r >> 17;
  // x_m_1 = 0x9; // = 2^3 + 1
  t = 0;
  asm("pclmulqdq %2, %1, %0;"
    : "=x"(t)
    : "x"(r_n_deg_r), "i"(0), "0"(x_m_1)
    );
  r = r_0_n_1 ^ t;
  t = r >> 17;
  r = r ^ (t << 3) ^ t;
  c_shares_p[1][0] = r & 0x1FFFF;

}

void CircuitContainer::squSharedExperimental33Verify(word* a_shares, word* b_shares, word* c_shares) {

  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;

  // Procedure for each party separately:
  // 1. Apply linear squaring algorithm
  // 2. Apply steps for modular reduction

  doubleword r;
  word u_0, u_1, u_2, u_3, c0, c1, T;
  // Party 1
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[0][0] & 0xFF;
  u_1 = (a_shares_p[0][0] & 0xFF00) >> 8;
  u_2 = (a_shares_p[0][0] & 0xFF0000) >> 16;
  u_3 = (a_shares_p[0][0] & 0xFF000000) >> 24;
  r = (this->squaring_precomp_[u_3] << 16) | this->squaring_precomp_[u_2];
  r = (r << 32) | (this->squaring_precomp_[u_1] << 16) | this->squaring_precomp_[u_0];
  // Last 1 bit (MS 33th bit)
  u_0 = (a_shares_p[0][0] & 0x100000000) >> 32;
  r = r | ((doubleword)(this->squaring_precomp_[u_0]) << 64);
  // Reduction
  c0 = r & 0x1FFFFFFFF; // LS 33 bits
  c1 = r >> 33; // MS 32 bits
  T = c1;
  c1 = c1 ^ (T >> 27) ^ (T >> 30); // 27 = 33 - 6, 30 = 33 - 3, 32 = 33 - 1 (omitted, all zeros), x^0 does not affect c1
  T = c1;
  c_shares_p[0][0] = c0 ^ ((T << 6) & 0x1FFFFFFFF) ^ ((T << 3) & 0x1FFFFFFFF) ^ ((T << 1) & 0x1FFFFFFFF) ^ T; // = c0, for x^6, x^3, x^1, x^0

  // Party 2
  // Linear Squaring
  // First 32 bits
  u_0 = a_shares_p[1][0] & 0xFF;
  u_1 = (a_shares_p[1][0] & 0xFF00) >> 8;
  u_2 = (a_shares_p[1][0] & 0xFF0000) >> 16;
  u_3 = (a_shares_p[1][0] & 0xFF000000) >> 24;
  r = (this->squaring_precomp_[u_3] << 16) | this->squaring_precomp_[u_2];
  r = (r << 32) | (this->squaring_precomp_[u_1] << 16) | this->squaring_precomp_[u_0];
  // Last 1 bit (MS 33th bit)
  u_0 = (a_shares_p[1][0] & 0x100000000) >> 32;
  r = r | ((doubleword)(this->squaring_precomp_[u_0]) << 64);
  // Reduction
  c0 = r & 0x1FFFFFFFF; // LS 33 bits
  c1 = r >> 33; // MS 32 bits
  T = c1;
  c1 = c1 ^ (T >> 27) ^ (T >> 30); // 27 = 33 - 6, 30 = 33 - 3, 32 = 33 - 1 (omitted, all zeros), x^0 does not affect c1
  T = c1;
  c_shares_p[1][0] = c0 ^ ((T << 6) & 0x1FFFFFFFF) ^ ((T << 3) & 0x1FFFFFFFF) ^ ((T << 1) & 0x1FFFFFFFF) ^ T; // = c0, for x^6, x^3, x^1, x^0

}

// word* shares contains the words of ALL shares (share_party_1_word_1, share_party_1_word_2, ... ,share_party_n_word_m-1, share_party_n_word_m) for n parties and m words
void CircuitContainer::mulSharedVerify(word* a_shares, word* b_shares, word* c_shares) {
  // How to:
  // 1. Compute value c_shares[0] using a_shares[0], a_shares[1], b_shares[0] and b_shares[1] (normally like above)
  // 2. Store value in View to c_shares[1] (this is the output that can't be computed here, so it has to be taken from the View)
  // All in all: c_shares[0] is computed like always, c_shares[1] gets the View values
  // Compute c_shares[0]
  word temp[this->gate_num_words_];
  uint32 offset = this->gate_num_words_;
  uint32 mul_gate_random_pointer = this->current_mul_gate_ * this->gate_num_words_;
  uint32 mul_gate_pointer_uchar = this->current_mul_gate_ * this->gate_size_;
  BigIntLib::Mul(temp, a_shares, b_shares);
  BigIntLib::Mul(c_shares, a_shares + offset, b_shares);
  BigIntLib::Add(temp, temp, c_shares);
  BigIntLib::Mul(c_shares, a_shares, b_shares + offset);
  BigIntLib::Add(temp, temp, c_shares);
  BigIntLib::Add(temp, temp, this->random_numbers_[0] + mul_gate_random_pointer);
  BigIntLib::Sub(c_shares, temp, this->random_numbers_[1] + mul_gate_random_pointer);
  // Get value from view and store in c_shares[1]
  memcpy(c_shares + offset, this->proof_->zs_[this->iteration_]->view_ + mul_gate_pointer_uchar, BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  // Write c_shares[0] to VerifyData (this is part of the View computed now)
  memcpy((this->verify_data_)->view_ + mul_gate_pointer_uchar, c_shares, BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  this->current_mul_gate_++;
}

void CircuitContainer::cubeSharedSign(word* a_shares, word* c_shares) {

  /*

  CBA = 000 (0)
  S(0, 0, 0) = [0 + 0 * 0, 0 + 0 * (0 + 0), 0 + 0 + 0 * (0 + 0)] = [0, 0, 0] -> should be [0, 0, 0]

  CBA = 001 (1)
  S(0, 0, 1) = [0 + 0 * 1, 0 + 1 * (0 + 0), 1 + 0 + 0 * (0 + 0)] = [0, 0, 1] -> should be [0, 0, 1]

  CBA = 010 (2)
  S(0, 1, 0) = [0 + 1 * 0, 1 + 0 * (1 + 0), 0 + 1 + 0 * (1 + 0)] = [0, 1, 1] -> should be [0, 1, 1]

  CBA = 011 (3)
  S(0, 1, 1) = [0 + 1 * 1, 1 + 1 * (1 + 0), 1 + 1 + 0 * (1 + 0)] = [1, 0, 0] -> should be [1, 0, 0]

  CBA = 100 (4)
  S(1, 0, 0) = [1 + 0 * 0, 0 + 0 * (0 + 1), 0 + 0 + 1 * (0 + 1)] = [1, 0, 1] -> should be [1, 0, 1]

  CBA = 101 (5)
  S(1, 0, 1) = [1 + 0 * 1, 0 + 1 * (0 + 1), 1 + 0 + 1 * (0 + 1)] = [1, 1, 0] -> should be [1, 1, 0]

  CBA = 110 (6)
  S(1, 1, 0) = [1 + 1 * 0, 1 + 0 * (1 + 1), 0 + 1 + 1 * (1 + 1)] = [1, 1, 1] -> should be [1, 1, 1]

  CBA = 111 (7)
  S(1, 1, 1) = [1 + 1 * 1, 1 + 1 * (1 + 1), 1 + 1 + 1 * (0 + 0)] = [0, 1, 0] -> should be [0, 1, 0]

  */

  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;

  uint32 mul_gate_random_pointer = this->current_mul_gate_ * this->gate_num_words_;
  uint32 mul_gate_pointer_uchar = this->current_mul_gate_ * this->gate_size_;

  // All values
  word c_0 = (a_shares_p[0][0] & 0x4) >> 2;
  word c_1 = (a_shares_p[1][0] & 0x4) >> 2;
  word c_2 = (a_shares_p[2][0] & 0x4) >> 2;
  word b_0 = (a_shares_p[0][0] & 0x2) >> 1;
  word b_1 = (a_shares_p[1][0] & 0x2) >> 1;
  word b_2 = (a_shares_p[2][0] & 0x2) >> 1;
  word a_0 = a_shares_p[0][0] & 0x1;
  word a_1 = a_shares_p[1][0] & 0x1;
  word a_2 = a_shares_p[2][0] & 0x1;

  // Random numbers
  word r_0_0 = (*(word*)(this->random_numbers_[0] + mul_gate_random_pointer) & 0x4) >> 2;
  word r_0_1 = (*(word*)(this->random_numbers_[0] + mul_gate_random_pointer) & 0x2) >> 1;
  word r_0_2 = *(word*)(this->random_numbers_[0] + mul_gate_random_pointer) & 0x1;
  word r_1_0 = (*(word*)(this->random_numbers_[1] + mul_gate_random_pointer) & 0x4) >> 2;
  word r_1_1 = (*(word*)(this->random_numbers_[1] + mul_gate_random_pointer) & 0x2) >> 1;
  word r_1_2 = *(word*)(this->random_numbers_[1] + mul_gate_random_pointer) & 0x1;
  word r_2_0 = (*(word*)(this->random_numbers_[2] + mul_gate_random_pointer) & 0x4) >> 2;
  word r_2_1 = (*(word*)(this->random_numbers_[2] + mul_gate_random_pointer) & 0x2) >> 1;
  word r_2_2 = *(word*)(this->random_numbers_[2] + mul_gate_random_pointer) & 0x1;


  c_shares_p[0][0] = 0;
  c_shares_p[1][0] = 0;
  c_shares_p[2][0] = 0;
  word t_mul = 0;

  word t_1 = (b_0 ^ c_0);
  word t_2 = (b_1 ^ c_1);
  word t_3 = (b_2 ^ c_2);

  // Unoptimized (same calculations performed multiple times, no De Morgan)
  // Party 0, C value
  t_mul = (b_0 & a_0) ^ (b_1 & a_0) ^ (b_0 & a_1) ^ r_0_0 ^ r_1_0;
  c_shares_p[0][0] = (c_0 ^ t_mul) << 2;

  // Party 0, B value
  t_mul = (a_0 & t_1) ^ (a_1 & t_1) ^ (a_0 & t_2) ^ r_0_1 ^ r_1_1;
  c_shares_p[0][0] = c_shares_p[0][0] | ((b_0 ^ t_mul) << 1);

  // Party 0, A value
  t_mul = (c_0 & t_1) ^ (c_1 & t_1) ^ (c_0 & t_2) ^ r_0_2 ^ r_1_2;
  c_shares_p[0][0] = c_shares_p[0][0] | (a_0 ^ b_0 ^ t_mul);

  // Copy shares for party 0
  memcpy((this->sign_data_)->views_[0] + mul_gate_pointer_uchar, c_shares_p[0], BigIntLib::field_size_bytes_);

  // Party 1, C value
  t_mul = (b_1 & a_1) ^ (b_2 & a_1) ^ (b_1 & a_2) ^ r_1_0 ^ r_2_0;
  c_shares_p[1][0] = (c_1 ^ t_mul) << 2;

  // Party 1, B value
  t_mul = (a_1 & t_2) ^ (a_2 & t_2) ^ (a_1 & t_3) ^ r_1_1 ^ r_2_1;
  c_shares_p[1][0] = c_shares_p[1][0] | ((b_1 ^ t_mul) << 1);

  // Party 1, A value
  t_mul = (c_1 & t_2) ^ (c_2 & t_2) ^ (c_1 & t_3) ^ r_1_2 ^ r_2_2;
  c_shares_p[1][0] = c_shares_p[1][0] | (a_1 ^ b_1 ^ t_mul);

  // Copy shares for party 1
  memcpy((this->sign_data_)->views_[1] + mul_gate_pointer_uchar, c_shares_p[1], BigIntLib::field_size_bytes_);

  // Party 2, C value
  t_mul = (b_2 & a_2) ^ (b_0 & a_2) ^ (b_2 & a_0) ^ r_2_0 ^ r_0_0;
  c_shares_p[2][0] = (c_2 ^ t_mul) << 2;

  // Party 2, B value
  t_mul = (a_2 & t_3) ^ (a_0 & t_3) ^ (a_2 & t_1) ^ r_2_1 ^ r_0_1;
  c_shares_p[2][0] = c_shares_p[2][0] | ((b_2 ^ t_mul) << 1);

  // Party 2, A value
  t_mul = (c_2 & t_3) ^ (c_0 & t_3) ^ (c_2 & t_1) ^ r_2_2 ^ r_0_2;
  c_shares_p[2][0] = c_shares_p[2][0] | (a_2 ^ b_2 ^ t_mul);

  // Copy shares for party 1
  memcpy((this->sign_data_)->views_[2] + mul_gate_pointer_uchar, c_shares_p[2], BigIntLib::field_size_bytes_);

  this->current_mul_gate_++;
}

void CircuitContainer::cubeSharedVerify(word* a_shares, word* c_shares) {

  /*
  S(C, B, A) = [C + B * (1 + A + C), C + A * (C + B), A + B + C * A]
  */

  word (*a_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) a_shares;
  word (*c_shares_p)[this->gate_num_words_] = (word (*)[this->gate_num_words_]) c_shares;

  uint32 offset = this->gate_num_words_;
  uint32 mul_gate_random_pointer = this->current_mul_gate_ * this->gate_num_words_;
  uint32 mul_gate_pointer_uchar = this->current_mul_gate_ * this->gate_size_;

  word c_0 = (a_shares_p[0][0] & 0x4) >> 2;
  word c_1 = (a_shares_p[1][0] & 0x4) >> 2;
  word b_0 = (a_shares_p[0][0] & 0x2) >> 1;
  word b_1 = (a_shares_p[1][0] & 0x2) >> 1;
  word a_0 = a_shares_p[0][0] & 0x1;
  word a_1 = a_shares_p[1][0] & 0x1;

  // Random numbers
  word r_0_0 = (*(word*)(this->random_numbers_[0] + mul_gate_random_pointer) & 0x4) >> 2;
  word r_0_1 = (*(word*)(this->random_numbers_[0] + mul_gate_random_pointer) & 0x2) >> 1;
  word r_0_2 = *(word*)(this->random_numbers_[0] + mul_gate_random_pointer) & 0x1;
  word r_1_0 = (*(word*)(this->random_numbers_[1] + mul_gate_random_pointer) & 0x4) >> 2;
  word r_1_1 = (*(word*)(this->random_numbers_[1] + mul_gate_random_pointer) & 0x2) >> 1;
  word r_1_2 = *(word*)(this->random_numbers_[1] + mul_gate_random_pointer) & 0x1;

  c_shares_p[0][0] = 0;
  c_shares_p[1][0] = 0;
  word t_mul = 0;

  word t_1 = (b_0 ^ c_0);

  // Without random values for testing, unoptimized (same calculations performed multiple times)
  // Party 0, C value
  t_mul = (b_0 & a_0) ^ (b_1 & a_0) ^ (b_0 & a_1) ^ r_0_0 ^ r_1_0;
  c_shares_p[0][0] = (c_0 ^ t_mul) << 2;

  // Party 0, B value
  t_mul = (a_0 & t_1) ^ (a_1 & t_1) ^ (a_0 & (c_1 ^ b_1)) ^ r_0_1 ^ r_1_1;
  c_shares_p[0][0] = c_shares_p[0][0] | ((b_0 ^ t_mul) << 1);

  // Party 0, A value
  t_mul = (c_0 & t_1) ^ (c_1 & t_1) ^ (c_0 & (b_1 ^ c_1)) ^ r_0_2 ^ r_1_2;
  c_shares_p[0][0] = (c_shares_p[0][0]) | (a_0 ^ b_0 ^ t_mul);

  // Get value from view and store in c_shares[1]
  memcpy(c_shares + offset, this->proof_->zs_[this->iteration_]->view_ + mul_gate_pointer_uchar, BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  // Write c_shares[0] to VerifyData (this is part of the View computed now)
  memcpy((this->verify_data_)->view_ + mul_gate_pointer_uchar, c_shares, BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_

  this->current_mul_gate_++;
}

// word* shares contains the words of ALL shares (share_party_1_word_1, share_party_1_word_2, ... ,share_party_n_word_m-1, share_party_n_word_m) for n parties and m words
void CircuitContainer::copyShares(word* from_shares, word* to_shares, uint32 party_size) {
  uint32 offset;
  for(uint32 i = 0; i < party_size; i++) {
    offset = i * this->gate_num_words_;
    memcpy(to_shares + offset, from_shares + offset, BigIntLib::field_size_bytes_); // Maybe use this->branch_size_ instead of this->gate_size_
  }
}

void CircuitContainer::prepareSharesFieldSign(uchar* x, word* value_shares_f, word* key_shares_f) {
  word (*value_shares_f_tmp)[this->party_size_][this->gate_num_words_] = (word (*)[this->party_size_][this->gate_num_words_]) value_shares_f;
  word (*key_shares_f_tmp)[this->party_size_][this->gate_num_words_] = (word (*)[this->party_size_][this->gate_num_words_]) key_shares_f;
  uint32 offset_1;
  word temp[this->gate_num_words_];
  memset(temp, 0, this->gate_size_);

  for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
    // Calculate last value in uchar and store all values into field values
    offset_1 = i * this->branch_size_;

    // Value shares (do 1 and 2 outside?)
    memcpy(value_shares_f_tmp[i][0], this->value_shares_[0] + offset_1, this->branch_size_);
    memset(value_shares_f_tmp[i][1], 0x0, this->branch_size_);
    memset(value_shares_f_tmp[i][2], 0x0, this->branch_size_);

    // Key shares
    memcpy(key_shares_f_tmp[i][0], this->key_shares_[0] + offset_1, this->branch_size_);
    memcpy(key_shares_f_tmp[i][1], this->key_shares_[1] + offset_1, this->branch_size_);
    memcpy(temp, this->key_ + offset_1, this->branch_size_);
    BigIntLib::TryReduce(key_shares_f_tmp[i][0]);
    BigIntLib::TryReduce(key_shares_f_tmp[i][1]);
    BigIntLib::TryReduce(temp);
    BigIntLib::Sub(key_shares_f_tmp[i][2], temp, key_shares_f_tmp[i][0]);
    BigIntLib::Sub(key_shares_f_tmp[i][2], key_shares_f_tmp[i][2], key_shares_f_tmp[i][1]);
    memcpy(this->key_shares_[2] + offset_1, key_shares_f_tmp[i][2], this->branch_size_);
  }
}

void CircuitContainer::prepareSharesFieldSignGeneric(uchar* x, word* value_shares_f, word* key_shares_f) {
  word (*value_shares_f_tmp)[this->party_size_][this->gate_num_words_] = (word (*)[this->party_size_][this->gate_num_words_]) value_shares_f;
  word (*key_shares_f_tmp)[this->party_size_][this->gate_num_words_] = (word (*)[this->party_size_][this->gate_num_words_]) key_shares_f;
  word temp_1[this->gate_num_words_];
  memset(temp_1, 0, this->gate_size_);

  uint32 words_full_count = ceil(this->value_size_ / 8.0);
  word* x_words = (word*)x;
  word temp_key[this->party_size_][words_full_count];
  memset(temp_key, 0, this->party_size_ * words_full_count * 8);
  memcpy(temp_key[0], this->key_shares_[0], this->key_size_);
  memcpy(temp_key[1], this->key_shares_[1], this->key_size_);
  memcpy(temp_key[2], this->key_, this->key_size_);
  word temp_2[words_full_count];
  memset(temp_2, 0, words_full_count * 8);

  uint32 num_bits;
  uint32 source_offset_words = 0;
  uint32 dest_offset_words = 0;
  uint32 dest_offset_bits = 0;
  word mask;
  uint32 bits_used = 0;
  uint32 bits_written;
  for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
    num_bits = BigIntLib::field_size_bits_;
    for(uint32 j = 0; j < this->gate_num_words_; j++) {
      bits_written = std::min(num_bits, (uint32)WORD_SIZE);
      mask = (0xFFFFFFFFFFFFFFFF >> (WORD_SIZE - bits_written));
      // Value
      // Share 1
      value_shares_f_tmp[i][0][j] = ((*(x_words + source_offset_words)) >> bits_used) & mask;

      // Share 2
      value_shares_f_tmp[i][1][j] = 0;

      // Share 3
      value_shares_f_tmp[i][2][j] = 0;

      // Key
      // Share 1
      key_shares_f_tmp[i][0][j] = ((*(temp_key[0] + source_offset_words)) >> bits_used) & mask;

      // Share 2
      key_shares_f_tmp[i][1][j] = ((*(temp_key[1] + source_offset_words)) >> bits_used) & mask;

      // Key k
      temp_1[j] = ((*(temp_key[2] + source_offset_words)) >> bits_used) & mask;

      bits_used += bits_written;

      if(bits_used >= WORD_SIZE) {
        bits_used -= WORD_SIZE;
        source_offset_words += 1;

        // Assign remaining bits_used bits from next word
        if(bits_used > 0) {
          // Valgrind warnings like "Invalid read of size 8" for the last word are to be expected here (e.g. MS 2 bytes are used as a whole word) (-> work with larger temp value instead and copy bytes at the end)
          value_shares_f_tmp[i][0][j] = (value_shares_f_tmp[i][0][j]) | ((((*(x_words + source_offset_words)) << (WORD_SIZE - bits_used)) >> (WORD_SIZE - bits_used)) << (bits_written - bits_used));
          key_shares_f_tmp[i][0][j] = (key_shares_f_tmp[i][0][j]) | ((((*(temp_key[0] + source_offset_words)) << (WORD_SIZE - bits_used)) >> (WORD_SIZE - bits_used)) << (bits_written - bits_used));
          key_shares_f_tmp[i][1][j] = (key_shares_f_tmp[i][1][j]) | ((((*(temp_key[1] + source_offset_words)) << (WORD_SIZE - bits_used)) >> (WORD_SIZE - bits_used)) << (bits_written - bits_used));
          temp_1[j] = (temp_1[j]) | ((((*(temp_key[2] + source_offset_words)) << (WORD_SIZE - bits_used)) >> (WORD_SIZE - bits_used)) << (bits_written - bits_used));
        }
      }

      num_bits -= bits_written;
    }

    BigIntLib::TryReduce(value_shares_f_tmp[i][0]);

    // Calc last share for key
    BigIntLib::TryReduce(key_shares_f_tmp[i][0]);
    BigIntLib::TryReduce(key_shares_f_tmp[i][1]);
    BigIntLib::TryReduce(temp_1);
    BigIntLib::Sub(key_shares_f_tmp[i][2], temp_1, key_shares_f_tmp[i][0]);
    BigIntLib::Sub(key_shares_f_tmp[i][2], key_shares_f_tmp[i][2], key_shares_f_tmp[i][1]);

    // Write last share to key_shares_tmp[2]
    num_bits = BigIntLib::field_size_bits_;
    word* pnt_temp_last_key = (word*)(temp_key[2]);
    for(uint32 j = 0; j < this->gate_num_words_; j++) {
      bits_written = std::min(num_bits, (uint32)WORD_SIZE);
      temp_2[dest_offset_words] = temp_2[dest_offset_words] | (key_shares_f_tmp[i][2][j] << dest_offset_bits); // no mask necessary, should be all zeros here
      num_bits -= bits_written;
      dest_offset_bits += bits_written;

      if(dest_offset_bits >= WORD_SIZE) {
        dest_offset_bits -= WORD_SIZE;
        dest_offset_words += 1;
        // Write remaining dest_offset_bits to next word
        uint32 shift = std::min(BigIntLib::field_size_bits_, (uint32)WORD_SIZE);
        if(dest_offset_bits > 0) {
          // Valgrind warnings like "Invalid read of size 8" for the last word are to be expected here (e.g. MS 2 bytes are used as a whole word) (-> work with larger temp value instead and copy bytes at the end)
          temp_2[dest_offset_words] = temp_2[dest_offset_words] | (key_shares_f_tmp[i][2][j] >> (shift - dest_offset_bits));
        }
      }
    }
  }

  memcpy(this->key_shares_[2], temp_2, this->key_size_);
}

void CircuitContainer::prepareSharesFieldVerify(word* value_shares_f, word* key_shares_f) {
  word (*value_shares_f_tmp)[this->party_size_ - 1][this->gate_num_words_] = (word (*)[this->party_size_ - 1][this->gate_num_words_]) value_shares_f;
  word (*key_shares_f_tmp)[this->party_size_ - 1][this->gate_num_words_] = (word (*)[this->party_size_ - 1][this->gate_num_words_]) key_shares_f;
  uint32 offset_1;
  for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
    // Calculate last value in uchar and store all values into field values
    offset_1 = i * this->branch_size_;

    // Value shares
    memcpy(value_shares_f_tmp[i][0], this->value_shares_[0] + offset_1, this->branch_size_);
    memcpy(value_shares_f_tmp[i][1], this->value_shares_[1] + offset_1, this->branch_size_);
    BigIntLib::TryReduce(value_shares_f_tmp[i][0]);
    BigIntLib::TryReduce(value_shares_f_tmp[i][1]);

    // Key shares
    memcpy(key_shares_f_tmp[i][0], this->key_shares_[0] + offset_1, this->branch_size_);
    memcpy(key_shares_f_tmp[i][1], this->key_shares_[1] + offset_1, this->branch_size_);
    BigIntLib::TryReduce(key_shares_f_tmp[i][0]);
    BigIntLib::TryReduce(key_shares_f_tmp[i][1]);
  }
}

void CircuitContainer::prepareSharesFieldVerifyGeneric(word* value_shares_f, word* key_shares_f) {
  word (*value_shares_f_tmp)[this->party_size_ - 1][this->gate_num_words_] = (word (*)[this->party_size_ - 1][this->gate_num_words_]) value_shares_f;
  word (*key_shares_f_tmp)[this->party_size_ - 1][this->gate_num_words_] = (word (*)[this->party_size_ - 1][this->gate_num_words_]) key_shares_f;

  uint32 words_full_count = ceil(this->value_size_ / 8.0);
  word temp_val[this->party_size_ - 1][words_full_count];
  word temp_key[this->party_size_ - 1][words_full_count];
  memset(temp_val, 0, (this->party_size_ - 1) * words_full_count * 8);
  memset(temp_key, 0, (this->party_size_ - 1) * words_full_count * 8);
  memcpy(temp_val[0], this->value_shares_[0], this->value_size_);
  memcpy(temp_val[1], this->value_shares_[1], this->value_size_);
  memcpy(temp_key[0], this->key_shares_[0], this->key_size_);
  memcpy(temp_key[1], this->key_shares_[1], this->key_size_);

  uint32 num_bits;
  uint32 source_offset_words = 0;
  word mask;
  uint32 bits_used = 0;
  uint32 bits_written;
  for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
    num_bits = BigIntLib::field_size_bits_;
    for(uint32 j = 0; j < this->gate_num_words_; j++) {
      bits_written = std::min(num_bits, (uint32)WORD_SIZE);
      mask = (0xFFFFFFFFFFFFFFFF >> (WORD_SIZE - bits_written));
      // Value
      // Share 1
      value_shares_f_tmp[i][0][j] = ((*(temp_val[0] + source_offset_words)) >> bits_used) & mask;

      // Share 2
      value_shares_f_tmp[i][1][j] = ((*(temp_val[1] + source_offset_words)) >> bits_used) & mask;

      // Key
      // Share 1
      key_shares_f_tmp[i][0][j] = ((*(temp_key[0] + source_offset_words)) >> bits_used) & mask;

      // Share 2
      key_shares_f_tmp[i][1][j] = ((*(temp_key[1] + source_offset_words)) >> bits_used) & mask;

      bits_used += bits_written;

      if(bits_used >= WORD_SIZE) {
        bits_used -= WORD_SIZE;
        source_offset_words += 1;

        // Assign remaining bits_used bits from next word
        if(bits_used > 0) {
          // Valgrind warnings like "Invalid read of size 8" for the last word are to be expected here (e.g. MS 2 bytes are used as a whole word) (-> work with larger temp value instead and copy bytes at the end)
          value_shares_f_tmp[i][0][j] = (value_shares_f_tmp[i][0][j]) | ((((*(temp_val[0] + source_offset_words)) << (WORD_SIZE - bits_used)) >> (WORD_SIZE - bits_used)) << (bits_written - bits_used));
          value_shares_f_tmp[i][1][j] = (value_shares_f_tmp[i][1][j]) | ((((*(temp_val[1] + source_offset_words)) << (WORD_SIZE - bits_used)) >> (WORD_SIZE - bits_used)) << (bits_written - bits_used));
          key_shares_f_tmp[i][0][j] = (key_shares_f_tmp[i][0][j]) | ((((*(temp_key[0] + source_offset_words)) << (WORD_SIZE - bits_used)) >> (WORD_SIZE - bits_used)) << (bits_written - bits_used));
          key_shares_f_tmp[i][1][j] = (key_shares_f_tmp[i][1][j]) | ((((*(temp_key[1] + source_offset_words)) << (WORD_SIZE - bits_used)) >> (WORD_SIZE - bits_used)) << (bits_written - bits_used));
        }
      }

      num_bits -= bits_written;
    }

    BigIntLib::TryReduce(value_shares_f_tmp[i][0]);
    BigIntLib::TryReduce(value_shares_f_tmp[i][1]);
    BigIntLib::TryReduce(key_shares_f_tmp[i][0]);
    BigIntLib::TryReduce(key_shares_f_tmp[i][1]);

  }
}

void CircuitContainer::outputSharesToBytes(word* output_shares, uint32 party_size) {
  word (*pointer)[party_size][this->gate_num_words_] = (word (*)[party_size][this->gate_num_words_]) output_shares;
  for(uint32 i = 0; i < party_size; i++) {
    for(uint32 j = 0; j < this->num_feistel_branches_; j++) {
      memcpy(this->value_shares_[i] + (j * this->branch_size_), pointer[(this->feistel_branch_indices_)[j]][i], this->branch_size_);
      //memcpy((this->sign_data_)->y_shares_ + (i * this->value_size_) + (j * this->branch_size_), pointer[(this->feistel_branch_indices_)[j]][i], this->branch_size_);
      //memcpy(this->value_shares_[i] + (j * this->branch_size_), output_shares + ((this->feistel_branch_indices_)[j]) * this->gate_num_words_ * party_size + i, this->branch_size_);

      //std::cout << "V1: " << std::hex << *(word*)(output_shares + ((this->feistel_branch_indices_)[j]) * this->gate_num_words_ * party_size + i) << std::dec << std::endl;
      //std::cout << "V2: " << std::hex << *(word*)(pointer[(this->feistel_branch_indices_)[j]][i]) << std::dec << std::endl;
    }
  }
}

void CircuitContainer::outputSharesToBytesGeneric(word* output_shares, uint32 party_size) {
  word (*pointer)[party_size][this->gate_num_words_] = (word (*)[party_size][this->gate_num_words_]) output_shares;

  uint32 words_full_count = ceil(this->value_size_ / 8.0);
  word temp[party_size][words_full_count];
  memset(temp, 0, party_size * words_full_count * 8);

  uint32 num_bits;
  uint32 dest_offset_words = 0;
  uint32 dest_offset_bits = 0;
  uint32 bits_written;
  for(uint32 i = 0; i < party_size; i++) {
    dest_offset_words = 0;
    dest_offset_bits = 0;
    for(uint32 j = 0; j < this->num_feistel_branches_; j++) {
      num_bits = BigIntLib::field_size_bits_;
      for(uint32 k = 0; k < this->gate_num_words_; k++) {
        bits_written = std::min(num_bits, (uint32)WORD_SIZE);
        *(temp[i] + dest_offset_words) = *(temp[i] + dest_offset_words) | (pointer[(this->feistel_branch_indices_)[j]][i][k] << dest_offset_bits); // no mask necessary, should be all zeros here
        num_bits -= bits_written;
        dest_offset_bits += bits_written;

        if(dest_offset_bits >= WORD_SIZE) {
          dest_offset_bits -= WORD_SIZE;
          dest_offset_words += 1;
          // Write remaining dest_offset_bits to next word
          uint32 shift = std::min(BigIntLib::field_size_bits_, (uint32)WORD_SIZE);
          if(dest_offset_bits > 0) {
            // Valgrind warnings like "Invalid read of size 8" for the last word are to be expected here (e.g. MS 2 bytes are used as a whole word) (-> work with larger temp value instead and copy bytes at the end)
            *(temp[i] + dest_offset_words) = *(temp[i] + dest_offset_words) | (pointer[(this->feistel_branch_indices_)[j]][i][k] >> (shift - dest_offset_bits));
          }
        }
      }
    }
    memcpy(this->value_shares_[i], temp[i], this->value_size_);
  }
}

void CircuitContainer::verifyCalcLastShare(uchar* y, word* value_shares_f) {
  word (*pointer)[this->party_size_ - 1][this->gate_num_words_] = (word (*)[this->party_size_ - 1][this->gate_num_words_]) value_shares_f;
  word words_temp[this->num_feistel_branches_][this->gate_num_words_];
  memset(words_temp, 0, this->num_feistel_branches_ * this->gate_size_);
  word value_temp[this->num_feistel_branches_][this->gate_num_words_];
  memset(value_temp, 0, this->num_feistel_branches_ * this->gate_size_);
  word given_temp[this->num_feistel_branches_][this->gate_num_words_];
  memset(given_temp, 0, this->num_feistel_branches_ * this->gate_size_);
  uint32 offset_1;
  for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
    offset_1 = i * this->branch_size_;
    memcpy(value_temp[i], y + offset_1, this->branch_size_);
    memcpy(given_temp[i], ((this->proof_)->zs_[this->iteration_])->y_share_ + offset_1, this->branch_size_);
    BigIntLib::Sub(words_temp[i], value_temp[i], pointer[(this->feistel_branch_indices_)[i]][0]); // First party is always the calculated party in the verify step
    BigIntLib::Sub(words_temp[i], words_temp[i], given_temp[i]);
    memcpy((this->verify_data_)->y_e2_ + offset_1, words_temp[i], this->branch_size_);
  }
}

void CircuitContainer::verifyCalcLastShareGeneric(uchar* y, word* value_shares_f) {
  word (*pointer)[this->party_size_ - 1][this->gate_num_words_] = (word (*)[this->party_size_ - 1][this->gate_num_words_]) value_shares_f;
  word words_temp[this->num_feistel_branches_][this->gate_num_words_];
  memset(words_temp, 0, this->num_feistel_branches_ * this->gate_size_);
  word value_temp[this->num_feistel_branches_][this->gate_num_words_];
  memset(value_temp, 0, this->num_feistel_branches_ * this->gate_size_);
  word given_temp[this->num_feistel_branches_][this->gate_num_words_];
  memset(given_temp, 0, this->num_feistel_branches_ * this->gate_size_);

  uint32 words_full_count = ceil(this->value_size_ / 8.0);
  word temp[3][words_full_count];
  memset(temp, 0, 3 * words_full_count * 8);
  memcpy(temp[0], y, this->value_size_);
  memcpy(temp[1], ((this->proof_)->zs_[this->iteration_])->y_share_, this->value_size_);

  uint32 num_bits;
  uint32 source_offset_words = 0;
  uint32 dest_offset_words = 0;
  uint32 dest_offset_bits = 0;
  word mask;
  uint32 bits_used = 0;
  uint32 bits_written;
  for(uint32 i = 0; i < this->num_feistel_branches_; i++) {
    num_bits = BigIntLib::field_size_bits_;
    for(uint32 j = 0; j < this->gate_num_words_; j++) {
      bits_written = std::min(num_bits, (uint32)WORD_SIZE);
      mask = (0xFFFFFFFFFFFFFFFF >> (WORD_SIZE - bits_written));
      // Value temp
      value_temp[i][j] = ((*(temp[0] + source_offset_words)) >> bits_used) & mask;

      // Given temp
      given_temp[i][j] = ((*(temp[1] + source_offset_words)) >> bits_used) & mask;

      bits_used += bits_written;

      if(bits_used >= WORD_SIZE) {
        bits_used -= WORD_SIZE;
        source_offset_words += 1;

        // Assign remaining bits_used bits from next word
        if(bits_used > 0) {
          // Valgrind warnings like "Invalid read of size 8" for the last word are to be expected here (e.g. MS 2 bytes are used as a whole word) (-> work with larger temp value instead and copy bytes at the end)
          value_temp[i][j] = (value_temp[i][j]) | ((((*(temp[0] + source_offset_words)) << (WORD_SIZE - bits_used)) >> (WORD_SIZE - bits_used)) << (bits_written - bits_used));
          given_temp[i][j] = (given_temp[i][j]) | ((((*(temp[1] + source_offset_words)) << (WORD_SIZE - bits_used)) >> (WORD_SIZE - bits_used)) << (bits_written - bits_used));
        }
      }

      num_bits -= bits_written;
    }

    // Calc last shares for value
    BigIntLib::Sub(words_temp[i], value_temp[i], pointer[(this->feistel_branch_indices_)[i]][0]); // First party is always the calculated party in the verify step
    BigIntLib::Sub(words_temp[i], words_temp[i], given_temp[i]);

    // Write last shares to (this->verify_data_)->y_e2_
    num_bits = BigIntLib::field_size_bits_;
    for(uint32 j = 0; j < this->gate_num_words_; j++) {
      bits_written = std::min(num_bits, (uint32)WORD_SIZE);
      //mask = (0xFFFFFFFFFFFFFFFF >> (WORD_SIZE - bits_written));
      *(temp[2] + dest_offset_words) = *(temp[2] + dest_offset_words) | (words_temp[i][j] << dest_offset_bits); // no mask necessary, should be all zeros here
      num_bits -= bits_written;
      dest_offset_bits += bits_written;

      if(dest_offset_bits >= WORD_SIZE) {
        dest_offset_bits -= WORD_SIZE;
        dest_offset_words += 1;
        // Write remaining dest_offset_bits to next word
        uint32 shift = std::min(BigIntLib::field_size_bits_, (uint32)WORD_SIZE);
        if(dest_offset_bits > 0) {
          // Valgrind warnings like "Invalid read of size 8" for the last word are to be expected here (e.g. MS 2 bytes are used as a whole word) (-> work with larger temp value instead and copy bytes at the end)
          *(temp[2] + dest_offset_words) = *(temp[2] + dest_offset_words) | (words_temp[i][j] >> (shift - dest_offset_bits));
        }
      }
    }
  }
  memcpy((this->verify_data_)->y_e2_, temp[2], this->value_size_);
}

uint64 CircuitContainer::rdtsc() {
  uint32 high, low;
  __asm__ __volatile__ ("rdtsc" : "=a" (low), "=d" (high));
  return (((uint64)high << 32) | low);
}
