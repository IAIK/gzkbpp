// Functions for operations in prime fields and binary fields
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

#include "BigIntLib.h"
#include <cmath> // log2(.)
#include <chrono>
#include <iostream>
#include <sstream>
#include <algorithm> // std::max, std::min
#include <iomanip> // For HEX output
#include <cstring> // memcpy
#include <openssl/rand.h> // For RAND_bytes(.)
#include <openssl/evp.h> // For AES

//#define VERBOSE
//#define VERBOSE_ERROR

uint32 BigIntLib::field_type_;
uint32 BigIntLib::field_size_bytes_;
uint32 BigIntLib::field_size_bits_;
uint32 BigIntLib::field_num_words_;
uint64 BigIntLib::msb_word_mask_;
uint32 BigIntLib::fastreduc_num_words_;
word* BigIntLib::modulo_;
EVP_CIPHER_CTX* BigIntLib::evp_cipher_ctx_;
uchar* BigIntLib::evp_cipher_buffer_;
int** BigIntLib::reduction_matrix_;
int** BigIntLib::reduction_matrix_positive_;
int** BigIntLib::reduction_matrix_negative_;
uint32 BigIntLib::mod_addition_weight_;
int** BigIntLib::mod_addition_matrix_;
uint32 BigIntLib::mod_subtraction_weight_;
int** BigIntLib::mod_subtraction_matrix_;

// Function variables
void (*ADD_FUNCTION)(word*, word*, word*);
void (*ADD_SPEC_FUNCTION)(word*, word*, word*);
void (*SUB_FUNCTION)(word*, word*, word*);
void (*MUL_FUNCTION)(word*, word*, word*);
void (*MUL_SPEC_FUNCTION)(word*, word*);
void (*REDUC_FUNCTION)(word*, word*);
void (*DOUBLE_ADD_FUNCTION)(word*, word*, word*);
void (*TIMES_2_FUNCTION)(word*, word*);
void (*TIMES_3_FUNCTION)(word*, word*);
bool (*GREATER_FUNCTION)(word*, word*);
void (*TRY_REDUCE_FUNCTION)(word*);
void (*GET_RANDOM_FIELD_ELEMENTS_FUNCTION)(void*, uchar*, uint32);
void (*GET_RANDOM_FIELD_ELEMENTS_NON_ZERO_FUNCTION)(void*, uchar*, uint32);

void BigIntLib::Init(uint32 branch_size_bits, uint32 field_type) {
  BigIntLib::field_type_ = field_type;
  BigIntLib::field_size_bits_ = branch_size_bits;

  // Size of the field in bits is branch_size_bits + 1 if it is a binary field
  /*
  if(field_type == 1) {
    BigIntLib::field_size_bits_ += 1;
  }
  */

  // Init buffer for cipher (AES in this case, 16 bytes output)
  BigIntLib::evp_cipher_buffer_ = new uchar[16];

  BigIntLib::field_size_bytes_ = ceil(((float)BigIntLib::field_size_bits_ / WORD_SIZE) * 8);
  BigIntLib::field_num_words_ = ceil((float)(BigIntLib::field_size_bytes_) / (WORD_SIZE / 8));

  uint32 msb_word_unused_bits = (WORD_SIZE - (BigIntLib::field_size_bits_ % WORD_SIZE)) % WORD_SIZE;
  BigIntLib::msb_word_mask_ = (0xFFFFFFFFFFFFFFFF >> msb_word_unused_bits);

  BigIntLib::modulo_ = new word[BigIntLib::field_num_words_];

  // This if construction is bad, change to something better later
  if(BigIntLib::field_type_ == 0) { // Prime field
    int* t_coeff;
    GREATER_FUNCTION = &BigIntLib::GreaterPF;
    TRY_REDUCE_FUNCTION = &BigIntLib::TryReducePF;
    GET_RANDOM_FIELD_ELEMENTS_FUNCTION = &BigIntLib::GetRandomFieldElementsPF;
    GET_RANDOM_FIELD_ELEMENTS_NON_ZERO_FUNCTION = &BigIntLib::GetRandomFieldElementsNonZeroPF;
    TIMES_2_FUNCTION = &BigIntLib::Times2PF;
    TIMES_3_FUNCTION = &BigIntLib::Times3PF;
    if(BigIntLib::field_size_bits_ == 3) {
      // Solinas prime 2^3 - 2 - 1
      BigIntLib::modulo_[0] = 0x5;
      // 2^3 - 2 - 1 = t^3 - t - 1, where t = 2^1
      BigIntLib::fastreduc_num_words_ = 3;
      t_coeff = new int[BigIntLib::fastreduc_num_words_]();
      t_coeff[0] = -1;
      t_coeff[1] = -1;

      ADD_FUNCTION = &BigIntLib::AddPF3;
      SUB_FUNCTION = &BigIntLib::SubPF3;
      MUL_FUNCTION = &BigIntLib::MulPF3Fast;
    }
    else if(BigIntLib::field_size_bits_ == 4) {
      // Solinas prime 2^4 - 2^2 - 1
      // Crandall prime 2^4 - 5 (same as Solinas)
      BigIntLib::modulo_[0] = 0xB;
      // 2^4 - 2^2 - 1 = t^2 - t - 1, where t = 2^2
      BigIntLib::fastreduc_num_words_ = 2;
      t_coeff = new int[BigIntLib::fastreduc_num_words_]();
      t_coeff[0] = -1;
      t_coeff[1] = -1;

      ADD_FUNCTION = &BigIntLib::AddPF4;
      SUB_FUNCTION = &BigIntLib::SubPF4;
      MUL_FUNCTION = &BigIntLib::MulPF4Fast;
      //MUL_FUNCTION = &BigIntLib::MulPF4FastCrandall;
      DOUBLE_ADD_FUNCTION = &BigIntLib::DoubleAddPF4;
    }
    else if(BigIntLib::field_size_bits_ == 16) {
      // Solinas prime 2^16 - 2^4 - 1
      // Crandall prime 2^32 - 17 (same as Solinas)
      BigIntLib::modulo_[0] = 0xFFEF;
      // 2^16 - 2^4 - 1 = t^4 - t - 1, where t = 2^4
      BigIntLib::fastreduc_num_words_ = 4;
      t_coeff = new int[BigIntLib::fastreduc_num_words_]();
      t_coeff[0] = -1;
      t_coeff[1] = -1;

      ADD_FUNCTION = &BigIntLib::AddPF16;
      SUB_FUNCTION = &BigIntLib::SubPF16;
      MUL_FUNCTION = &BigIntLib::MulPF16Fast;
      //MUL_FUNCTION = &BigIntLib::MulPF16FastCrandall;

      // Experimental
      ADD_SPEC_FUNCTION = &BigIntLib::AddSpecPF16;
      MUL_SPEC_FUNCTION = &BigIntLib::MulSpecPF16;
      REDUC_FUNCTION = &BigIntLib::SolinasReducPF16;
    }
    else if(BigIntLib::field_size_bits_ == 32) {
      // Solinas prime 2^32 - 2^4 - 1
      // Crandall prime 2^32 - 5
      BigIntLib::modulo_[0] = 0xFFFFFFEF; // Solinas USE THIS!
      //BigIntLib::modulo_[0] = 0xFEFFFFFF; // 2^32 - 2^24 - 1 = t^4 - t^3 - 1 ONLY FOR PRINTING OUT MATRICES!
      //BigIntLib::modulo_[0] = 0xFFFFFFFB; // Crandall
      // 2^32 - 2^4 - 1 = t^8 - t - 1, where t = 2^4
      BigIntLib::fastreduc_num_words_ = 8; // USE THIS!
      //BigIntLib::fastreduc_num_words_ = 4; // 2^32 - 2^24 - 1 = t^4 - t^3 - 1 ONLY FOR PRINTING OUT MATRICES!
      t_coeff = new int[BigIntLib::fastreduc_num_words_]();
      t_coeff[0] = -1;
      t_coeff[1] = -1; // USE THIS!
      //t_coeff[3] = -1; // 2^32 - 2^24 - 1 = t^4 - t^3 - 1 ONLY FOR PRINTING OUT MATRICES!


      ADD_FUNCTION = &BigIntLib::AddPF32;
      SUB_FUNCTION = &BigIntLib::SubPF32;
      MUL_FUNCTION = &BigIntLib::MulPF32Fast;
      //MUL_FUNCTION = &BigIntLib::MulPF32FastCrandall;
    }
    else if(BigIntLib::field_size_bits_ == 64) {
      // Solinas prime 2^64 - 2^8 - 1
      // Crandall prime 2^64 - 59
      BigIntLib::modulo_[0] = 0xFFFFFFFFFFFFFEFF; // Solinas
      //BigIntLib::modulo_[0] = 0xFFFFFFFFFFFFFFC5; // Crandall
      // 2^64 - 2^8 - 1 = t^8 - t - 1, where t = 2^8
      BigIntLib::fastreduc_num_words_ = 8;
      t_coeff = new int[BigIntLib::fastreduc_num_words_]();
      t_coeff[0] = -1;
      t_coeff[1] = -1;

      ADD_FUNCTION = &BigIntLib::AddPF64;
      SUB_FUNCTION = &BigIntLib::SubPF64;
      //MUL_FUNCTION = &BigIntLib::MulPF64;
      MUL_FUNCTION = &BigIntLib::MulPF64Fast;
      //MUL_FUNCTION = &BigIntLib::MulPF64FastCrandall;
    }
    else if(BigIntLib::field_size_bits_ == 136) {
      // Solinas prime 2^136 - 2^8 - 1 (alternative with 16-bit-words: 2^144 - 2^128 - 1)
      BigIntLib::modulo_[0] = 0xFFFFFFFFFFFFFEFF;
      BigIntLib::modulo_[1] = 0xFFFFFFFFFFFFFFFF;
      BigIntLib::modulo_[2] = 0x00000000000000FF;
      // 2^136 - 2^8 - 1 = t^17 - t - 1, where t = 2^8
      BigIntLib::fastreduc_num_words_ = 17;
      t_coeff = new int[BigIntLib::fastreduc_num_words_]();
      t_coeff[0] = -1;
      t_coeff[1] = -1;

      ADD_FUNCTION = &BigIntLib::AddPF136;
      SUB_FUNCTION = &BigIntLib::SubPF136;
      MUL_FUNCTION = &BigIntLib::MulPF136Fast;
    }
    else if(BigIntLib::field_size_bits_ == 256) {
      // Implement
      // 2^256 - 2^184 + 2^32 + 1 = t^32 - t^23 + t^4 + 1, where t = 2^8 (Generalized Mersenne Prime)
      BigIntLib::modulo_[0] = 0x0000000100000001;
      BigIntLib::modulo_[1] = 0x0000000000000000;
      BigIntLib::modulo_[2] = 0xFF00000000000000;
      BigIntLib::modulo_[3] = 0xFFFFFFFFFFFFFFFF;
      BigIntLib::fastreduc_num_words_ = 32;
      t_coeff = new int[BigIntLib::fastreduc_num_words_]();
      t_coeff[0] = 1;
      t_coeff[4] = 1;
      t_coeff[23] = -1;
      // TESTING with t^7 - t^3 + 1
      // BigIntLib::field_size_bytes_ = 7;
      //t_coeff[0] = 1;
      //t_coeff[3] = -1;
      // TESTING with t^16 - t - 1
      // BigIntLib::field_size_bytes_ = 16;
      //t_coeff[0] = -1;
      //t_coeff[1] = -1;
      // TESTING with t^14 - t^7 - 1
      //BigIntLib::field_size_bytes_ = 14;
      //t_coeff[0] = -1;
      //t_coeff[7] = -1;

      ADD_FUNCTION = &BigIntLib::AddPF256;
      SUB_FUNCTION = &BigIntLib::SubPF256;
      MUL_FUNCTION = &BigIntLib::MulPF256Fast;
    }
    else if(BigIntLib::field_size_bits_ == 272) {
      // 272-bit gates like MiMC+, write functions for this...
      BigIntLib::modulo_[0] = 0xFFFFFEFFFFFFFFFF;
      BigIntLib::modulo_[1] = 0xFFFFFFFFFFFFFFFF;
      BigIntLib::modulo_[2] = 0xFFFFFFFFFFFFFFFF;
      BigIntLib::modulo_[3] = 0xFFFFFFFFFFFFFFFF;
      BigIntLib::modulo_[4] = 0x000000000000FFFF;
      // 2^272 - 2^40 - 1 = t^34 - t^5 - 1, where t = 2^8 (this is a Solinas prima)
      BigIntLib::fastreduc_num_words_ = 34;
      t_coeff = new int[BigIntLib::fastreduc_num_words_]();
      t_coeff[0] = -1;
      t_coeff[5] = -1;

      ADD_FUNCTION = &BigIntLib::AddPF272;
      SUB_FUNCTION = &BigIntLib::SubPF272;
      MUL_FUNCTION = &BigIntLib::MulPF272Fast;
    }
    else {
      std::cout << "[BigIntLib] Error: Field size " << BigIntLib::field_size_bits_ << " not supported for this field type." << std::endl;
      exit(1);
    }

    // Compute matrices used for reduction
    BigIntLib::ComputeReductionMatrix(t_coeff);
    BigIntLib::ComputeModularAdditionMatrix();
    BigIntLib::ComputeModularSubtractionMatrix();

    // Clean up
    delete[] t_coeff;
  }
  else if(BigIntLib::field_type_ == 1) { // Binary field
    GREATER_FUNCTION = &BigIntLib::GreaterBF;
    TRY_REDUCE_FUNCTION = &BigIntLib::TryReduceBF;
    GET_RANDOM_FIELD_ELEMENTS_FUNCTION = &BigIntLib::GetRandomFieldElementsBF;
    GET_RANDOM_FIELD_ELEMENTS_NON_ZERO_FUNCTION = &BigIntLib::GetRandomFieldElementsNonZeroBF;
    if(BigIntLib::field_size_bits_ == 3) {
      // Irreducible polynomial x^3 + x + 1
      BigIntLib::modulo_[0] = 0xB;

      ADD_FUNCTION = &BigIntLib::XorBF3;
      SUB_FUNCTION = &BigIntLib::XorBF3;
      MUL_FUNCTION = &BigIntLib::MulBF3Fast;
      TIMES_2_FUNCTION = &BigIntLib::Times2BF3;
      TIMES_3_FUNCTION = &BigIntLib::Times3BF3;
    }
    else if(BigIntLib::field_size_bits_ == 17) {
      // Irreducible polynomial x^17 + x^3 + 1
      BigIntLib::modulo_[0] = 0x20009;

      ADD_FUNCTION = &BigIntLib::XorBF17;
      SUB_FUNCTION = &BigIntLib::XorBF17;
      MUL_FUNCTION = &BigIntLib::MulBF17Fast;
      TIMES_2_FUNCTION = &BigIntLib::Times2BF17;
      TIMES_3_FUNCTION = &BigIntLib::Times3BF17;
    }
    else if(BigIntLib::field_size_bits_ == 33) {
      // Irreducible polynomial x^33 + x^6 + x^3 + x + 1
      BigIntLib::modulo_[0] = 0x20000004B;

      ADD_FUNCTION = &BigIntLib::XorBF33;
      SUB_FUNCTION = &BigIntLib::XorBF33;
      MUL_FUNCTION = &BigIntLib::MulBF33Fast;
      TIMES_2_FUNCTION = &BigIntLib::Times2BF33;
      TIMES_3_FUNCTION = &BigIntLib::Times3BF33;
    }
    else if(BigIntLib::field_size_bits_ == 65) {
      // Irreducible polynomial x^65 + x^4 + x^3 + x + 1
      BigIntLib::modulo_[0] = 0x1B;
      BigIntLib::modulo_[1] = 0x2;

      ADD_FUNCTION = &BigIntLib::XorBF65;
      SUB_FUNCTION = &BigIntLib::XorBF65;
      MUL_FUNCTION = &BigIntLib::MulBF65Fast;
      TIMES_2_FUNCTION = &BigIntLib::Times2BF65;
      TIMES_3_FUNCTION = &BigIntLib::Times3BF65;
    }
    else {
      std::cout << "[BigIntLib] Error: Field size " << BigIntLib::field_size_bits_ << " not supported for this field type." << std::endl;
      exit(1);
    }
  }
  else {
    std::cout << "[BigIntLib] Error: Please choose field type 0 (prime field) or 1 (binary field)." << std::endl;
    exit(1);
  }
  
  // TODO TESTING
  //BigIntLib::field_size_bytes_ = 3;
  //memset(t_coeff, 0, sizeof(int) * field_size);
  //for(uint32 i = 0; i < 34; i++) t_coeff[i] = 0;
  //t_coeff[0] = 1;
  //t_coeff[1] = -1;

  // Init AES
  BigIntLib::evp_cipher_ctx_ = EVP_CIPHER_CTX_new();
}

void BigIntLib::CleanUp() {
  delete[] BigIntLib::modulo_;

  // Matrices (this can be shortened in a single for loop), do it only for prime fields
  if(BigIntLib::field_type_ == 0) {
    uint32 n = BigIntLib::fastreduc_num_words_;
    for(uint32 i = 0; i < n; i++) {
      delete[] BigIntLib::reduction_matrix_[i];
    }
    delete[] BigIntLib::reduction_matrix_;
    for(uint32 i = 0; i < n; i++) {
      delete[] BigIntLib::reduction_matrix_positive_[i];
    }
    delete[] BigIntLib::reduction_matrix_positive_;
    for(uint32 i = 0; i < n; i++) {
      delete[] BigIntLib::reduction_matrix_negative_[i];
    }
    delete[] BigIntLib::reduction_matrix_negative_;
    for(uint32 i = 0; i < BigIntLib::mod_addition_weight_; i++) {
      delete[] BigIntLib::mod_addition_matrix_[i];
    }
    delete[] BigIntLib::mod_addition_matrix_;
    for(uint32 i = 0; i < BigIntLib::mod_subtraction_weight_; i++) {
      delete[] BigIntLib::mod_subtraction_matrix_[i];
    }
    delete[] BigIntLib::mod_subtraction_matrix_;
  }

  // Clear AES and buffer
  delete[] BigIntLib::evp_cipher_buffer_;
  EVP_CIPHER_CTX_free(BigIntLib::evp_cipher_ctx_);
}

void BigIntLib::Add(word* c, word* a, word* b) {
  #ifdef VERBOSE
  if(BigIntLib::field_type_ == 0) {
    std::cout << "--- ADD ---" << std::endl;
    std::cout << "Copy for Sage:" << std::endl;
    std::cout << "(0x" << BigIntLib::ToString(a) << " + 0x" << BigIntLib::ToString(b) << ") \% 0x" << BigIntLib::ToString(BigIntLib::modulo_) << std::endl;
  }
  else if(BigIntLib::field_type_ == 1) {
    std::cout << "--- ADD ---" << std::endl;
    std::cout << "Copy for Sage:" << std::endl;
    std::cout << "0x" << BigIntLib::ToString(a) << " ^ 0x" << BigIntLib::ToString(b) << std::endl;
  }
  #endif
  ADD_FUNCTION(c, a, b);
  #ifdef VERBOSE
  std::cout << "Result: 0x" << BigIntLib::ToString(c) << std::endl;
  #endif

  #ifdef VERBOSE_ERROR
  if(BigIntLib::field_type_ == 0) {
    if(BigIntLib::Greater(a, BigIntLib::modulo_)) {
      std::cout << "[ADD] ERROR!!! a > modulo" << std::endl;
      exit(1);
    }

    if(BigIntLib::Greater(b, BigIntLib::modulo_)) {
      std::cout << "[ADD] ERROR!!! b > modulo" << std::endl;
      exit(1);
    }

    if(BigIntLib::Greater(c, BigIntLib::modulo_)) {
      std::cout << "[ADD] ERROR!!! Result > modulo" << std::endl;
      std::cout << "--- FAILED ADD ---" << std::endl;
      std::cout << "Copy for Sage:" << std::endl;
      std::cout << "(0x" << BigIntLib::ToString(a) << " + 0x" << BigIntLib::ToString(b) << ") \% 0x" << BigIntLib::ToString(BigIntLib::modulo_) << std::endl;
      std::cout << "Result: 0x" << BigIntLib::ToString(c) << std::endl;
      exit(1);
    }
  }
  #endif
}

void BigIntLib::AddSpec(word* c, word* a, word* b) {
  ADD_SPEC_FUNCTION(c, a, b);
}

void BigIntLib::Sub(word* c, word* a, word* b) {
  #ifdef VERBOSE
  if(BigIntLib::field_type_ == 0) {
    std::cout << "--- SUB ---" << std::endl;
    std::cout << "Copy for Sage:" << std::endl;
    std::cout << "(0x" << BigIntLib::ToString(a) << " - 0x" << BigIntLib::ToString(b) << ") \% 0x" << BigIntLib::ToString(BigIntLib::modulo_) << std::endl;
  }
  else if(BigIntLib::field_type_ == 1) {
    std::cout << "--- SUB ---" << std::endl;
    std::cout << "Copy for Sage:" << std::endl;
    std::cout << "0x" << BigIntLib::ToString(a) << " ^ 0x" << BigIntLib::ToString(b) << std::endl;
  }
  #endif
  SUB_FUNCTION(c, a, b);
  #ifdef VERBOSE
  std::cout << "Result: 0x" << BigIntLib::ToString(c) << std::endl;
  #endif

  #ifdef VERBOSE_ERROR
  if(BigIntLib::field_type_ == 0) {
    if(BigIntLib::Greater(a, BigIntLib::modulo_)) {
      std::cout << "[SUB] ERROR!!! a > modulo" << std::endl;
      exit(1);
    }

    if(BigIntLib::Greater(b, BigIntLib::modulo_)) {
      std::cout << "[SUB] ERROR!!! b > modulo" << std::endl;
      exit(1);
    }

    if(BigIntLib::Greater(c, BigIntLib::modulo_)) {
      std::cout << "[SUB] ERROR!!! Result > modulo" << std::endl;
      std::cout << "--- FAILED SUB ---" << std::endl;
      std::cout << "Copy for Sage:" << std::endl;
      std::cout << "(0x" << BigIntLib::ToString(a) << " - 0x" << BigIntLib::ToString(b) << ") \% 0x" << BigIntLib::ToString(BigIntLib::modulo_) << std::endl;
      std::cout << "Result: 0x" << BigIntLib::ToString(c) << std::endl;
      exit(1);
    }
  }
  #endif
  
}

void BigIntLib::Mul(word* c, word* a, word* b) {
  #ifdef VERBOSE
  if(BigIntLib::field_type_ == 0) {
    std::cout << "--- MUL ---" << std::endl;
    std::cout << "Copy for Sage:" << std::endl;
    std::cout << "(0x" << BigIntLib::ToString(a) << " * 0x" << BigIntLib::ToString(b) << ") \% 0x" << BigIntLib::ToString(BigIntLib::modulo_) << std::endl;
  }
  else if(BigIntLib::field_type_ == 1) {
    std::cout << "--- MUL ---" << std::endl;
    std::cout << "Copy for Sage:" << std::endl;
    std::cout << "FIELD_SIZE = " << BigIntLib::field_size_bits_ << std::endl;
    std::cout << "a = 0x" << BigIntLib::ToString(a) << std::endl;
    std::cout << "b = 0x" << BigIntLib::ToString(b) << std::endl;
  }
  #endif
  MUL_FUNCTION(c, a, b);
  #ifdef VERBOSE
  std::cout << "Result: 0x" << BigIntLib::ToString(c) << std::endl;
  #endif
}

void BigIntLib::MulSpec(word* c, word* a) {
  MUL_SPEC_FUNCTION(c, a);
}

void BigIntLib::SolinasReduc(word* out, word* in) {
  REDUC_FUNCTION(out, in);
}

void BigIntLib::DoubleAdd(word* c, word* a, word* b) {
  #ifdef VERBOSE
  std::cout << "--- DOUBLEADD ---" << std::endl;
  std::cout << "Copy for Sage:" << std::endl;
  std::cout << "(0x" << BigIntLib::ToString(a) << " * 0x" << BigIntLib::ToString(b) << ") \% 0x" << BigIntLib::ToString(BigIntLib::modulo_) << std::endl;
  #endif
  DOUBLE_ADD_FUNCTION(c, a, b);
  #ifdef VERBOSE
  std::cout << "Result: 0x" << BigIntLib::ToString(c) << std::endl;
  #endif
}

void BigIntLib::Times2(word* c, word* a) {
  TIMES_2_FUNCTION(c, a);
}

void BigIntLib::Times3(word* c, word* a) {
  TIMES_3_FUNCTION(c, a);
}

void BigIntLib::FillRandom(void* destination, uchar* seed, uint32 num_bytes) {
  // This method assumes num_bytes mod 16 = 0
  // 128 bit IV
  uint32 aes_blocksize = 16;
  static const unsigned char iv[16] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '1', '2', '3', '4', '5'};
  EVP_EncryptInit_ex(BigIntLib::evp_cipher_ctx_, EVP_aes_128_ctr(), NULL, seed, iv);

  static const unsigned char plaintext[16] = {'0'};

  /*
  uint32 num_executions = num_bytes / aes_blocksize;
  int length;
  uchar* dest = (uchar*)destination;
  for(uint32 i = 0; i < num_executions; i++) {
    EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, dest, &length, plaintext, sizeof(plaintext));
    dest += aes_blocksize;
  }
  */

  uint32 count = num_bytes;
  uchar* dest = (uchar*)destination;
  int length;
  for(; count >= aes_blocksize; count -= aes_blocksize, dest += aes_blocksize) {
    EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, dest, &length, plaintext, sizeof(plaintext));
  }

  if(count) {
    EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, dest, &length, plaintext, count);
  }

  //memset(destination, 0, num_bytes); // FOR TESTING, TEMP
}

void BigIntLib::Print(void* source, uint32 num_bytes) {
  // TEMP
  uchar* pointer = (uchar*)source;
  for(uint32 i = 0; i < num_bytes; i++) {
    std::cout << std::setfill('0') << std::setw(2) << std::hex << (uint32)(pointer[num_bytes - i - 1]);
  }
  std::cout << std::dec << std::endl;
}

std::string BigIntLib::ToString(void* source, uint32 num_bytes) {
  if(num_bytes == 0) num_bytes = BigIntLib::field_size_bytes_;
  uchar* pointer = (uchar*)source;
  std::ostringstream string_stream;
  for(uint32 i = 0; i < num_bytes; i++) {
    string_stream << std::setfill('0') << std::setw(2) << std::hex << (uint32)(pointer[num_bytes - i - 1]);
  }
  std::string ret_string = string_stream.str();
  return ret_string;
}

void BigIntLib::ComputeReductionMatrix(int* t_coeff) {
  // Method due to Jerome A. Solinas, "Generalized Mersenne Numbers", 1999
  // The result is a n x n matrix of int numbers, where n is the number of terms (depends on word size)
  // E.g.: 2^64 - 2^8 - 1 = t^8 - t - 1, where t = 2^8
  uint32 n = BigIntLib::fastreduc_num_words_;
  //int t_coeff[n]; // TODO Hard-coded
  //memset(t_coeff, 0, sizeof(int) * n); // TODO Hard-coded
  //t_coeff[0] = -1; // TODO Hard-coded
  //t_coeff[1] = -1; // TODO Hard-coded

  // Allocate memory
  BigIntLib::reduction_matrix_ = new int*[n];
  for(uint32 i = 0; i < n; i++) {
    BigIntLib::reduction_matrix_[i] = new int[n];
  }
  BigIntLib::reduction_matrix_positive_ = new int*[n];
  for(uint32 i = 0; i < n; i++) {
    BigIntLib::reduction_matrix_positive_[i] = new int[n];
  }
  BigIntLib::reduction_matrix_negative_ = new int*[n];
  for(uint32 i = 0; i < n; i++) {
    BigIntLib::reduction_matrix_negative_[i] = new int[n];
  }

  for(uint32 j = 0; j < n; j++) {
    BigIntLib::reduction_matrix_[0][j] = -(t_coeff[j]);
  }
  for(uint32 i = 1; i < n; i++) {
    for(uint32 j = 0; j < n; j++) {
      if(j == 0) {
        BigIntLib::reduction_matrix_[i][j] = -(t_coeff[0]) * BigIntLib::reduction_matrix_[i - 1][n - 1];
      }
      else {
        BigIntLib::reduction_matrix_[i][j] = BigIntLib::reduction_matrix_[i - 1][j - 1] - (t_coeff[j] * BigIntLib::reduction_matrix_[i - 1][n - 1]);
      }
    }
  }

  // Split in positive and negative part
  int val_temp;
  for(uint32 i = 0; i < n; i++) {
    for(uint32 j = 0; j < n; j++) {
      val_temp = BigIntLib::reduction_matrix_[i][j];
      if(val_temp > 0) {
        BigIntLib::reduction_matrix_positive_[i][j] = val_temp;
        BigIntLib::reduction_matrix_negative_[i][j] = 0;
      }
      else if(val_temp < 0) {
        BigIntLib::reduction_matrix_positive_[i][j] = 0;
        BigIntLib::reduction_matrix_negative_[i][j] = val_temp;
      }
      else {
        BigIntLib::reduction_matrix_positive_[i][j] = 0;
        BigIntLib::reduction_matrix_negative_[i][j] = 0;
      }
      //if(BigIntLib::reduction_matrix_positive_[i][j] == 2) BigIntLib::reduction_matrix_positive_[i][j] = 1; // TEMP TESTING
    }
  }

  // Get modular addition weight and modular subtraction weight
  int max_addition = 0;
  int current_addition;
  int max_subtraction = 0;
  int current_subtraction;
  for(uint32 j = 0; j < n; j++) {
    current_addition = 0;
    current_subtraction = 0;
    for(uint32 i = 0; i < n; i++) {
      current_addition += BigIntLib::reduction_matrix_positive_[i][j];
      current_subtraction += BigIntLib::reduction_matrix_negative_[i][j];
    }
    if(current_addition > max_addition) max_addition = current_addition;
    if(current_subtraction < max_subtraction) max_subtraction = current_subtraction;
  }
  BigIntLib::mod_addition_weight_ = max_addition;
  BigIntLib::mod_subtraction_weight_ = max_subtraction * (-1);

  #ifdef VERBOSE
  // Print reduction matrix
  std::cout << "Computed reduction matrix X:" << std::endl;
  for(uint32 i = 0; i < n; i++) {
    for(uint32 j = 0; j < n; j++) {
       std::cout << BigIntLib::reduction_matrix_[i][j] << "  ";
    }
    std::cout << std::endl;
  }

  // Print positive part
  std::cout << "Positive part X+ of reduction matrix X:" << std::endl;
  for(uint32 i = 0; i < n; i++) {
    for(uint32 j = 0; j < n; j++) {
       std::cout << BigIntLib::reduction_matrix_positive_[i][j] << "  ";
    }
    std::cout << std::endl;
  }

  // Print Negative part
  std::cout << "Negative part X- of reduction matrix X:" << std::endl;
  for(uint32 i = 0; i < n; i++) {
    for(uint32 j = 0; j < n; j++) {
       std::cout << BigIntLib::reduction_matrix_negative_[i][j] << "  ";
    }
    std::cout << std::endl;
  }

  std::cout << "Modular addition weight: " << BigIntLib::mod_addition_weight_ << std::endl;
  std::cout << "Modular subtraction weight: " << BigIntLib::mod_subtraction_weight_ << std::endl;
  #endif

}

void BigIntLib::ComputeModularAdditionMatrix() {
  // Method due to Jerome A. Solinas, "Generalized Mersenne Numbers", 1999
  // Specific algorithm due to Mario Taschwer, "Modular Multiplication Using Special Prime Moduli", ‎2001
  //uint32 mod_addition_weight = 3; // TODO Hard-coded
  uint32 mod_weight = BigIntLib::mod_addition_weight_;
  uint32 n = BigIntLib::fastreduc_num_words_;

  // Allocate memory
  BigIntLib::mod_addition_matrix_ = new int*[mod_weight];
  for(uint32 i = 0; i < mod_weight; i++) {
    BigIntLib::mod_addition_matrix_[i] = new int[n];
  }
  
  int i;
  int val_temp;
  for(uint32 k = 0; k < mod_weight; k++) {
    for(uint32 j = 0; j < n; j++) {
      i = 0;
      while(i < n && BigIntLib::reduction_matrix_positive_[i][j] == 0) {
        i++;
      }
      if(i < n) {
        val_temp = n + i;
        if(val_temp == (n * 2)) val_temp = -1; // TEMP
        BigIntLib::mod_addition_matrix_[k][j] = val_temp;
        BigIntLib::reduction_matrix_positive_[i][j]--;
      }
      else {
        val_temp = 2 * n;
        if(val_temp == (n * 2)) val_temp = -1; // TEMP
        BigIntLib::mod_addition_matrix_[k][j] = val_temp;
      }
    }
  }

  #ifdef VERBOSE
  // Print modular addition matrix
  std::cout << "Computed modular addition matrix A+:" << std::endl;
  for(uint32 i = 0; i < mod_weight; i++) {
    for(uint32 j = 0; j < n; j++) {
      std::cout << BigIntLib::mod_addition_matrix_[i][j] << "  ";
    }
    std::cout << std::endl;
  }
  #endif
}

void BigIntLib::ComputeModularSubtractionMatrix() {
  // Method due to Jerome A. Solinas, "Generalized Mersenne Numbers", 1999
  //uint32 mod_addition_weight = 3; // TODO Hard-coded
  uint32 mod_weight = BigIntLib::mod_subtraction_weight_;
  uint32 n = BigIntLib::fastreduc_num_words_;

  // Allocate memory
  BigIntLib::mod_subtraction_matrix_ = new int*[mod_weight];
  for(uint32 i = 0; i < mod_weight; i++) {
    BigIntLib::mod_subtraction_matrix_[i] = new int[n];
  }
  
  int i;
  int val_temp;
  for(uint32 k = 0; k < mod_weight; k++) {
    for(uint32 j = 0; j < n; j++) {
      i = 0;
      while(i < n && BigIntLib::reduction_matrix_negative_[i][j] == 0) {
        i++;
      }
      if(i < n) {
        val_temp = n + i;
        if(val_temp == (n * 2)) val_temp = -1; // TEMP
        BigIntLib::mod_subtraction_matrix_[k][j] = val_temp;
        BigIntLib::reduction_matrix_negative_[i][j]++; // Difference to the calculation of the modular addition matrix
      }
      else {
        val_temp = 2 * n;
        if(val_temp == (n * 2)) val_temp = -1; // TEMP
        BigIntLib::mod_subtraction_matrix_[k][j] = val_temp;
      }
    }
  }

  #ifdef VERBOSE
  // Print modular addition matrix
  std::cout << "Computed modular subtraction matrix A-:" << std::endl;
  for(uint32 i = 0; i < mod_weight; i++) {
    for(uint32 j = 0; j < n; j++) {
       std::cout << BigIntLib::mod_subtraction_matrix_[i][j] << "  ";
    }
    std::cout << std::endl;
  }
  #endif
}

bool BigIntLib::IsZero(word* a) {
  // memcmp with test block (all zeros)?
  for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
    if(a[i] != 0) return false;
  }
  return true;
}

bool BigIntLib::Smaller(word* a, word* b) {
  for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
    if(a[BigIntLib::field_num_words_ - i - 1] < b[BigIntLib::field_num_words_ - i - 1]) return true;
    else if(a[BigIntLib::field_num_words_ - i - 1] > b[BigIntLib::field_num_words_ - i - 1]) return false;
  }
  return false;
}

bool BigIntLib::Greater(word* a, word* b) { // Change and rename to GreaterEqual later!
  return GREATER_FUNCTION(a, b);
}

bool BigIntLib::GreaterPF(word* a, word* b) { // Change and rename to GreaterEqualPF later!
  for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
    if(a[BigIntLib::field_num_words_ - i - 1] > b[BigIntLib::field_num_words_ - i - 1]) return true;
    else if(a[BigIntLib::field_num_words_ - i - 1] < b[BigIntLib::field_num_words_ - i - 1]) return false;
  }
  return false;
}

bool BigIntLib::GreaterBF(word* a, word* b) { // Change and rename to GreaterEqualBF later!
  return false; // TODO
}

uchar BigIntLib::Compare(word* a, word* b) {
  // 0 .. a < b
  // 1 .. a = b
  // 2 .. a > b
  for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
    if(a[BigIntLib::field_num_words_ - i - 1] < b[BigIntLib::field_num_words_ - i - 1]) return 0;
    else if(a[BigIntLib::field_num_words_ - i - 1] > b[BigIntLib::field_num_words_ - i - 1]) return 2;
  }
  return 1;
}

void BigIntLib::TryReduce(word* a) { // Checks if value is in field, reduces if not
  TRY_REDUCE_FUNCTION(a);
}

void BigIntLib::TryReducePF(word* a) { // Checks if value is in field, reduces if not
  if(BigIntLib::GreaterPF(a, BigIntLib::modulo_)) BigIntLib::Sub(a, a, BigIntLib::modulo_);
}

void BigIntLib::TryReduceBF(word* a) { // Checks if value is in field, reduces if not TODO
  return;
}

void BigIntLib::GetRandomFieldElements(void* destination, uchar* seed, uint32 num_elements) {
  GET_RANDOM_FIELD_ELEMENTS_FUNCTION(destination, seed, num_elements);
}

void BigIntLib::GetRandomFieldElementsNonZero(void* destination, uchar* seed, uint32 num_elements) {
  GET_RANDOM_FIELD_ELEMENTS_NON_ZERO_FUNCTION(destination, seed, num_elements);
}

void BigIntLib::GetRandomFieldElementsPF(void* destination, uchar* seed, uint32 num_elements) {
  uchar* dest_pointer = (uchar*)destination;

  // Init AES
  uint32 aes_blocksize = 16;
  int length;
  static const unsigned char iv[16] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '1', '2', '3', '4', '5'};
  EVP_EncryptInit_ex(BigIntLib::evp_cipher_ctx_, EVP_aes_128_ctr(), NULL, seed, iv);
  static const unsigned char plaintext[16] = {'0'};

  EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
  uint32 buffer_bytes_left = 16;
  uint32 buffer_bytes_used = 0;
  bool current_element_in_field = false;
  uint32 gate_size = BigIntLib::field_num_words_ * 8;
  uint32 msb_word_index = BigIntLib::field_num_words_ - 1;
  //uint32 bytes_for_element = BigIntLib::field_size_bytes_; // rounded up to the next larger byte value for BigIntLib::field_size_bits_ mod 8 != 0
  uint32 bytes_still_to_write;
  uint32 bytes_write_now;
  //uint32 written;

  uchar* current_write_pointer_base;
  uchar* current_write_pointer;
  word* current_element_number;

  // For each element
  for(uint32 i = 0; i < num_elements; i++) {
    current_write_pointer_base = dest_pointer + i * gate_size;
    current_element_number = (word*)current_write_pointer_base;
    //current_write_pointer = dest_pointer + i * gate_size;
    // This writes one complete element
    do {
      bytes_still_to_write = BigIntLib::field_size_bytes_;
      current_write_pointer = current_write_pointer_base;
      //current_element_number = (word*)current_write_pointer;
      //written = 0;
      //while(written < bytes_for_element) { // better: while(bytes_still_to_write > 0) and remove the "written" variable
      while(bytes_still_to_write > 0) {
        /*
        buffer_bytes_used = 16 - buffer_bytes_left;
        if(bytes_still_to_write <= buffer_bytes_left) {
          // There is enough bytes left
          memcpy(current_write_pointer, BigIntLib::evp_cipher_buffer_ + buffer_bytes_used, bytes_still_to_write);
          buffer_bytes_left -= bytes_still_to_write;
          bytes_still_to_write = 0;
          //written += bytes_still_to_write;
        }
        else if(buffer_bytes_left > 0) {
          // There is not enough bytes left, but still some amount left (write it)
          memcpy(current_write_pointer, BigIntLib::evp_cipher_buffer_ + buffer_bytes_used, buffer_bytes_left);
          //written += buffer_bytes_left;
          bytes_still_to_write -= buffer_bytes_left;
          current_write_pointer += buffer_bytes_left;
          // Get new AES bytes
          EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
          buffer_bytes_left = 16;
        }
        else {
          // Buffer is empty, get new AES bytes (maybe write next value already here?)
          EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
          buffer_bytes_left = 16;
        }
        */

        //std::cout << "bytes_still_to_write: " << bytes_still_to_write << std::endl;
        //std::cout << "buffer_bytes_left: " << buffer_bytes_left << std::endl;
        //std::cout << "buffer_bytes_used: " << buffer_bytes_used << std::endl;
        bytes_write_now = std::min(bytes_still_to_write, buffer_bytes_left);
        memcpy(current_write_pointer, BigIntLib::evp_cipher_buffer_ + buffer_bytes_used, bytes_write_now);
        buffer_bytes_left -= bytes_write_now;
        bytes_still_to_write -= bytes_write_now;
        current_write_pointer += bytes_write_now;
        buffer_bytes_used += bytes_write_now;
        if(buffer_bytes_left == 0) {
          // Buffer is empty, get new AES bytes (maybe write next value already here?)
          EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
          buffer_bytes_left = 16;
          buffer_bytes_used = 0;
        }

        /*
        std::cout << "buffer_bytes_left: " << buffer_bytes_left << std::endl;
        std::cout << "buffer_bytes_used: " << buffer_bytes_used << std::endl;
        std::cout << "bytes_still_to_write: " << bytes_still_to_write << std::endl;
        */
      }

      //exit(1);

      // Cancel out last word
      current_element_number[msb_word_index] &= BigIntLib::msb_word_mask_;
      current_element_in_field = BigIntLib::Smaller(current_element_number, BigIntLib::modulo_);
      //std::cout << "msb_word_index: " << msb_word_index << std::endl;
      //BigIntLib::Print(current_element_number, gate_size);

      /*
      if(current_element_in_field == false) {
        std::cout << "current_write_pointer: " << std::hex << (uint64)current_write_pointer << std::dec << std::endl;
        std::cout << "NOT IN FIELD!!! Value:" << std::endl;
        BigIntLib::Print(current_element_number, gate_size);
        //exit(1);
      }
      */
      
    } while(!current_element_in_field);
  }
}

void BigIntLib::GetRandomFieldElementsBF(void* destination, uchar* seed, uint32 num_elements) {
  uchar* dest_pointer = (uchar*)destination;

  // Init AES
  uint32 aes_blocksize = 16;
  int length;
  static const unsigned char iv[16] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '1', '2', '3', '4', '5'};
  EVP_EncryptInit_ex(BigIntLib::evp_cipher_ctx_, EVP_aes_128_ctr(), NULL, seed, iv);
  static const unsigned char plaintext[16] = {'0'};

  EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
  uint32 buffer_bytes_left = 16;
  uint32 buffer_bytes_used = 0;
  uint32 gate_size = BigIntLib::field_num_words_ * 8;
  uint32 msb_word_index = BigIntLib::field_num_words_ - 1;
  //uint32 bytes_for_element = BigIntLib::field_size_bytes_; // rounded up to the next larger byte value for BigIntLib::field_size_bits_ mod 8 != 0
  uint32 bytes_still_to_write;
  uint32 bytes_write_now;

  uchar* current_write_pointer_base;
  uchar* current_write_pointer;
  word* current_element_number;

  for(uint32 i = 0; i < num_elements; i++) {
    current_write_pointer_base = dest_pointer + i * gate_size;
    current_element_number = (word*)current_write_pointer_base;
    bytes_still_to_write = BigIntLib::field_size_bytes_;
    current_write_pointer = current_write_pointer_base;
    // This writes one complete element
    while(bytes_still_to_write > 0) {
      /*
      if(bytes_still_to_write <= buffer_bytes_left) {
        // There is enough bytes left
        memcpy(current_write_pointer, BigIntLib::evp_cipher_buffer_ + buffer_bytes_used, bytes_still_to_write);
        buffer_bytes_left -= bytes_still_to_write;
        bytes_still_to_write = 0;
      }
      else if(buffer_bytes_left > 0) {
        // There is not enough bytes left, but still some amount left (write it)
        memcpy(current_write_pointer, BigIntLib::evp_cipher_buffer_ + buffer_bytes_used, buffer_bytes_left);
        bytes_still_to_write -= buffer_bytes_left;
        current_write_pointer += buffer_bytes_left;
        // Get new AES bytes
        EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
        buffer_bytes_left = 16;
      }
      else {
        // Buffer is empty, get new AES bytes
        EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
        buffer_bytes_left = 16;
      }
      */

      bytes_write_now = std::min(bytes_still_to_write, buffer_bytes_left);
      memcpy(current_write_pointer, BigIntLib::evp_cipher_buffer_ + buffer_bytes_used, bytes_write_now);
      buffer_bytes_left -= bytes_write_now;
      bytes_still_to_write -= bytes_write_now;
      current_write_pointer += bytes_write_now;
      buffer_bytes_used += bytes_write_now;
      if(buffer_bytes_left == 0) {
        // Buffer is empty, get new AES bytes (maybe write next value already here?)
        EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
        buffer_bytes_left = 16;
        buffer_bytes_used = 0;
      }

    }

    // Cancel out last word
    current_element_number[msb_word_index] &= BigIntLib::msb_word_mask_;
  }
}

void BigIntLib::GetRandomFieldElementsNonZeroPF(void* destination, uchar* seed, uint32 num_elements) {
  uchar* dest_pointer = (uchar*)destination;

  // Init AES
  uint32 aes_blocksize = 16;
  int length;
  static const unsigned char iv[16] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '1', '2', '3', '4', '5'};
  EVP_EncryptInit_ex(BigIntLib::evp_cipher_ctx_, EVP_aes_128_ctr(), NULL, seed, iv);
  static const unsigned char plaintext[16] = {'0'};

  EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
  uint32 buffer_bytes_left = 16;
  uint32 buffer_bytes_used = 0;
  bool current_element_in_field = false;
  uint32 gate_size = BigIntLib::field_num_words_ * 8;
  uint32 msb_word_index = BigIntLib::field_num_words_ - 1;
  uint32 bytes_still_to_write;
  uint32 bytes_write_now;

  uchar* current_write_pointer_base;
  uchar* current_write_pointer;
  word* current_element_number;

  // For each element
  for(uint32 i = 0; i < num_elements; i++) {
    current_write_pointer_base = dest_pointer + i * gate_size;
    current_element_number = (word*)current_write_pointer_base;
    // This writes one complete element
    do {
      bytes_still_to_write = BigIntLib::field_size_bytes_;
      current_write_pointer = current_write_pointer_base;
      while(bytes_still_to_write > 0) {
        bytes_write_now = std::min(bytes_still_to_write, buffer_bytes_left);
        memcpy(current_write_pointer, BigIntLib::evp_cipher_buffer_ + buffer_bytes_used, bytes_write_now);
        buffer_bytes_left -= bytes_write_now;
        bytes_still_to_write -= bytes_write_now;
        current_write_pointer += bytes_write_now;
        buffer_bytes_used += bytes_write_now;

        if(buffer_bytes_left == 0) {
          // Buffer is empty, get new AES bytes (maybe write next value already here?)
          EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
          buffer_bytes_left = 16;
          buffer_bytes_used = 0;
        }

      }

      // Cancel out last word
      current_element_number[msb_word_index] &= BigIntLib::msb_word_mask_;
      current_element_in_field = BigIntLib::Smaller(current_element_number, BigIntLib::modulo_);
      
    } while(!current_element_in_field || BigIntLib::IsZero(current_element_number));
  }
}

void BigIntLib::GetRandomFieldElementsNonZeroBF(void* destination, uchar* seed, uint32 num_elements) {
  uchar* dest_pointer = (uchar*)destination;

  // Init AES
  uint32 aes_blocksize = 16;
  int length;
  static const unsigned char iv[16] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '1', '2', '3', '4', '5'};
  EVP_EncryptInit_ex(BigIntLib::evp_cipher_ctx_, EVP_aes_128_ctr(), NULL, seed, iv);
  static const unsigned char plaintext[16] = {'0'};

  EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
  uint32 buffer_bytes_left = 16;
  uint32 buffer_bytes_used = 0;
  uint32 gate_size = BigIntLib::field_num_words_ * 8;
  uint32 msb_word_index = BigIntLib::field_num_words_ - 1;
  //uint32 bytes_for_element = BigIntLib::field_size_bytes_; // rounded up to the next larger byte value for BigIntLib::field_size_bits_ mod 8 != 0
  uint32 bytes_still_to_write;
  uint32 bytes_write_now;

  uchar* current_write_pointer_base;
  uchar* current_write_pointer;
  word* current_element_number;

  for(uint32 i = 0; i < num_elements; i++) {
    current_write_pointer_base = dest_pointer + i * gate_size;
    current_element_number = (word*)current_write_pointer_base;
    // This writes one complete element
    do {
      bytes_still_to_write = BigIntLib::field_size_bytes_;
      current_write_pointer = current_write_pointer_base;
      while(bytes_still_to_write > 0) {
        bytes_write_now = std::min(bytes_still_to_write, buffer_bytes_left);
        memcpy(current_write_pointer, BigIntLib::evp_cipher_buffer_ + buffer_bytes_used, bytes_write_now);
        buffer_bytes_left -= bytes_write_now;
        bytes_still_to_write -= bytes_write_now;
        current_write_pointer += bytes_write_now;
        buffer_bytes_used += bytes_write_now;

        if(buffer_bytes_left == 0) {
          // Buffer is empty, get new AES bytes (maybe write next value already here?)
          EVP_EncryptUpdate(BigIntLib::evp_cipher_ctx_, BigIntLib::evp_cipher_buffer_, &length, plaintext, sizeof(plaintext));
          buffer_bytes_left = 16;
          buffer_bytes_used = 0;
        }

      }

      // Cancel out last word
      current_element_number[msb_word_index] &= BigIntLib::msb_word_mask_;
    } while(BigIntLib::IsZero(current_element_number));
  }
}

void BigIntLib::AddPF3(word* c, word* a, word* b) {
  word r = *a + *b;
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];
  *c = r;
}

void BigIntLib::SubPF3(word* c, word* a, word* b) {
  word r = (*a > *b) ? (*a - *b) : (*a + BigIntLib::modulo_[0] - *b);
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];
  *c = r;
}

void BigIntLib::MulPF3Fast(word* c, word* a, word* b) {
  uchar r = (*a * *b);

  // Modular additions
  // 3  3  4
  // 5  4  5
  // -1  5  -1

  uchar temp_var = 0x0;
  *c = r & 0x7;
  //temp_var = ((r & 0x8) >> 3) | ((r & 0x18) >> 2);
  *c = *c + (((r & 0x8) >> 3) | ((r & 0x18) >> 2));
  if(*c >= BigIntLib::modulo_[0]) *c -= BigIntLib::modulo_[0];
  //temp_var = ((r & 0x20) >> 5) | ((r & 0x30) >> 3);
  *c = *c + (((r & 0x20) >> 5) | ((r & 0x30) >> 3));
  if(*c >= BigIntLib::modulo_[0]) *c -= BigIntLib::modulo_[0];
  //temp_var = ((r & 0x20) >> 4);
  *c = *c + (((r & 0x20) >> 4));
  if(*c >= BigIntLib::modulo_[0]) *c -= BigIntLib::modulo_[0];
  
}

void BigIntLib::AddPF4(word* c, word* a, word* b) {
  word r = *a + *b;
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];
  *c = r;
}

void BigIntLib::SubPF4(word* c, word* a, word* b) {
  word r = (*a > *b) ? (*a - *b) : (*a + BigIntLib::modulo_[0] - *b);
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];
  *c = r;
}

void BigIntLib::MulPF4Fast(word* c, word* a, word* b) {
  uchar r = (*a * *b);

  // Modular additions
  // 2  2
  // 3  3
  // -1  3

  uchar temp_var = 0x0;
  *c = r & 0xF;
  temp_var = (((uchar)(r << 2) >> 6) << 2) | ((uchar)(r << 2) >> 6);
  *c = *c + temp_var;
  if(*c >= BigIntLib::modulo_[0]) *c -= BigIntLib::modulo_[0];
  temp_var = ((r >> 6) << 2) | (r >> 6);
  *c = *c + temp_var;
  if(*c >= BigIntLib::modulo_[0]) *c -= BigIntLib::modulo_[0];
  temp_var = (r >> 6) << 2;
  *c = *c + temp_var;
  if(*c >= BigIntLib::modulo_[0]) *c -= BigIntLib::modulo_[0];
}

void BigIntLib::MulPF4FastCrandall(word* c, word* a, word* b) {
  word r = (*a * *b);
  word c_small = 5;
  word q_i = r >> 4;
  r = r & 0xF;
  word t;

  // WHILE 1.1
  t = c_small * q_i;
  r = r + (t & 0xF) * (q_i > 0);
  q_i = t >> 4;

  // WHILE 1.2
  t = c_small * q_i;
  r = r + (t & 0xF) * (q_i > 0);

  // WHILE 2.1
  r -= BigIntLib::modulo_[0] * (r >= BigIntLib::modulo_[0]);

  // WHILE 2.2
  r -= BigIntLib::modulo_[0] * (r >= BigIntLib::modulo_[0]);
  
  *c = r;
}

void BigIntLib::DoubleAddPF4(word* c, word* a, word* b) {
  // Adds both 4-bit sides of the bytes separately
  // Least significant 4 bits
  word r1 = (*a & 0xF) + (*b & 0xF);
  if(r1 >= BigIntLib::modulo_[0]) r1 -= BigIntLib::modulo_[0];
  // Most significant 4 bits
  word r2 = ((*a & 0xF0) >> 4) + ((*b & 0xF0) >> 4);
  if(r2 >= BigIntLib::modulo_[0]) r2 -= BigIntLib::modulo_[0];
  // Assign
  *c = (r2 << 4) | r1;
}

void BigIntLib::AddPF16(word* c, word* a, word* b) {
  word r = *a + *b;
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];
  *c = r;
}

void BigIntLib::AddSpecPF16(word* c, word* a, word* b) {
  // No reduction
  *c = *a + *b;
}

void BigIntLib::SubPF16(word* c, word* a, word* b) {
  word r = (*a > *b) ? (*a - *b) : (*a + BigIntLib::modulo_[0] - *b);
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];
  *c = r;
}

void BigIntLib::MulPF16Fast(word* c, word* a, word* b) {
  word r = (*a * *b);

  // Modular additions
  // 4  4  5  6
  // 7  5  6  7
  // -1  7  -1  -1

  
  word temp_var;
  word t = r & 0xFFFF;
  temp_var = ((r & 0x000F0000) >> 16) | ((r & 0x0FFF0000) >> 12);
  t = t + temp_var;
  if(t >= BigIntLib::modulo_[0]) t -= BigIntLib::modulo_[0];
  temp_var = ((r & 0xF0000000) >> 28) | ((r & 0xFFF00000) >> 16);
  t = t + temp_var;
  if(t >= BigIntLib::modulo_[0]) t -= BigIntLib::modulo_[0];
  temp_var = ((r & 0xF0000000) >> 24);
  t = t + temp_var;
  if(t >= BigIntLib::modulo_[0]) t -= BigIntLib::modulo_[0];
  *c = t;
}

void BigIntLib::MulSpecPF16(word* c, word* a) {
  // Constant modular multiplication with 16
  word r = (*a << 4);

  // r is something like 0x000xxxxx (MS 12 bits are 0)

  word temp_var;
  word t = r & 0xFFFF;
  temp_var = ((r & 0x000F0000) >> 16) | ((r & 0x000F0000) >> 12);
  t = t + temp_var;
  if(t >= BigIntLib::modulo_[0]) t -= BigIntLib::modulo_[0];
  // Since the MS 12 bits of r are 0, the following steps of the Solinas reduction are not necessary
  *c = t;
}

void BigIntLib::SolinasReducPF16(word* out, word* in) {
  word r = *in;
  word temp_var;
  word t = r & 0xFFFF;
  temp_var = ((r & 0x000F0000) >> 16) | ((r & 0x0FFF0000) >> 12);
  t = t + temp_var;
  if(t >= BigIntLib::modulo_[0]) t -= BigIntLib::modulo_[0];
  temp_var = ((r & 0xF0000000) >> 28) | ((r & 0xFFF00000) >> 16);
  t = t + temp_var;
  if(t >= BigIntLib::modulo_[0]) t -= BigIntLib::modulo_[0];
  temp_var = ((r & 0xF0000000) >> 24);
  t = t + temp_var;
  if(t >= BigIntLib::modulo_[0]) t -= BigIntLib::modulo_[0];
  *out = t;
}

void BigIntLib::MulPF16FastCrandall(word* c, word* a, word* b) {
  word r = (*a * *b);
  word c_small = 17;
  word q_i = r >> 16;
  r = r & 0xFFFF;
  word t;

  // WHILE 1.1
  t = c_small * q_i;
  r = r + (t & 0xFFFF) * (q_i > 0);
  q_i = t >> 16;

  // WHILE 1.2
  t = c_small * q_i;
  r = r + (t & 0xFFFF) * (q_i > 0);

  // WHILE 2.1
  r -= BigIntLib::modulo_[0] * (r >= BigIntLib::modulo_[0]);

  // WHILE 2.2
  r -= BigIntLib::modulo_[0] * (r >= BigIntLib::modulo_[0]);

  *c = r;
}

void BigIntLib::AddPF32(word* c, word* a, word* b) {
  word r = (*a) + (*b);
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];
  *c = r;
}

void BigIntLib::SubPF32(word* c, word* a, word* b) {
  word r = (*a > *b) ? (*a - *b) : (*a + BigIntLib::modulo_[0] - *b);
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];
  *c = r;
}

void BigIntLib::MulPF32Fast(word* c, word* a, word* b) {
  word r = (*a * *b);

  // Modular additions
  // 8  8  9  10  11  12  13  14
  // 15  9  10  11  12  13  14  15
  // -1  15  -1  -1  -1  -1  -1  -1

  /*
  word temp_var;
  word x_tmp;
  x_tmp = r & 0xFFFFFFFF;
  
  temp_var = ((r >> 32) & 0xF) | (((r >> 32) & 0x0FFFFFFF) << 4);
  x_tmp = x_tmp + temp_var;
  if(x_tmp >= BigIntLib::modulo_[0]) x_tmp -= BigIntLib::modulo_[0];
  temp_var = ((r >> 60) & 0xF) | (((r >> 36) & 0x0FFFFFFF) << 4);
  x_tmp = x_tmp + temp_var;
  if(x_tmp >= BigIntLib::modulo_[0]) x_tmp -= BigIntLib::modulo_[0];
  temp_var = (((r >> 60) & 0xF) << 4);
  x_tmp = x_tmp + temp_var;
  if(x_tmp >= BigIntLib::modulo_[0]) x_tmp -= BigIntLib::modulo_[0];
  *c = (word)x_tmp;
  */

  // If computation is wrong, check first commented (old) lines!!!
  word temp_var;
  word t;
  t = r & 0xFFFFFFFF;
  //temp_var = ((r >> 32) & 0xF) | (((r >> 32) & 0x0FFFFFFF) << 4);
  temp_var = ((r >> 32) & 0xF) | ((r & 0xFFFFFFF00000000) >> 28);
  t = t + temp_var;
  if(t >= BigIntLib::modulo_[0]) t -= BigIntLib::modulo_[0];
  //t -= BigIntLib::modulo_[0] * (t >= BigIntLib::modulo_[0]);
  //temp_var = ((r >> 60) & 0xF) | (((r >> 36) & 0x0FFFFFFF) << 4);
  temp_var = (r >> 60) | ((r & 0xFFFFFFF000000000) >> 32);
  //temp_var = (r >> 60) | ((r >> 36) << 4);
  t = t + temp_var;
  if(t >= BigIntLib::modulo_[0]) t -= BigIntLib::modulo_[0];
  //t -= BigIntLib::modulo_[0] * (t >= BigIntLib::modulo_[0]);
  //temp_var = (((r >> 60) & 0xF) << 4);
  temp_var = ((r >> 60) << 4);
  t = t + temp_var;
  if(t >= BigIntLib::modulo_[0]) t -= BigIntLib::modulo_[0];
  //t -= BigIntLib::modulo_[0] * (t >= BigIntLib::modulo_[0]);
  *c = t;
}

void BigIntLib::MulPF32FastCrandall(word* c, word* a, word* b) {
  word r = (*a * *b);
  word c_small = 5; // 2^2 + 1 from p_32 = 2^32 - 2^2 - 1 = 2^32 - c
  word q_i = r >> 32;
  r = r & 0xFFFFFFFF;
  word t;

  // c * x = (c - 1) * x + x

  // WHILE 1.1
  //if(q_i > 0) { // This if-clause could be omitted with very high probability
    t = c_small * q_i;
    //t = (q_i << 2) + q_i;
    //r = r + (t & 0xFFFFFFFF);
    r = r + (t & 0xFFFFFFFF) * (q_i > 0);
    q_i = t >> 32;
  //}

  // WHILE 1.2
  //if(q_i > 0) { // This if-clause cannot be omitted
    t = c_small * q_i;
    //t = (q_i << 2) + q_i;
    //r = r + (t & 0xFFFFFFFF);
    r = r + (t & 0xFFFFFFFF) * (q_i > 0);
    //q_i = t >> 32;
  //}

  // WHILE 2.1
  /*
  if(r >= BigIntLib::modulo_[0]) // This if-clause cannot be omitted
    r -= BigIntLib::modulo_[0];
  // WHILE 2.2
  if(r >= BigIntLib::modulo_[0]) // This if-clause cannot be omitted
    r -= BigIntLib::modulo_[0];
  */
  // WHILE 2.1
  r -= BigIntLib::modulo_[0] * (r >= BigIntLib::modulo_[0]);

  // WHILE 2.2
  r -= BigIntLib::modulo_[0] * (r >= BigIntLib::modulo_[0]);
  
  *c = r;
}

void BigIntLib::AddPF64(word* c, word* a, word* b) {
  //*c = ((doubleword)(*a) + (doubleword)(*b)) % BigIntLib::modulo_[0];
  doubleword r = (doubleword)(*a) + *b;
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];
  //if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];
  *c = (word)r;
}

void BigIntLib::SubPF64(word* c, word* a, word* b) {
  //*c = (*a > *b) ? (*a - *b) : (*a + BigIntLib::modulo_[0] - *b);
  //if(*c >= BigIntLib::modulo_[0]) *c -= BigIntLib::modulo_[0];
  doubleword r = (*a > *b) ? (*a - *b) : (*a + BigIntLib::modulo_[0] - *b);
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];
  *c = (word)r;
}

void BigIntLib::MulPF64(word* c, word* a, word* b) {
  //std::cout << "--- MUL ---" << std::endl;
  //std::cout << "Copy for WolframAlpha:" << std::endl;
  //std::cout << "(" << *a << " * " << *b << ") mod " << BigIntLib::modulo_[0] << std::endl;
  *c = ((doubleword)(*a) * (doubleword)(*b)) % BigIntLib::modulo_[0];
  //std::cout << "Result: " << *c << std::endl;
}

void BigIntLib::MulPF64Fast(word* c, word* a, word* b) {
  // TODO Everything is hard-coded here, change it
  // Method due to Jerome A. Solinas, "Generalized Mersenne Numbers", 1999
  // Inputs are matrices X and A, which are called the reduction matrix and the modular addition matrix, respectively.
  // These are calculated in the beginning. The additions below can easily be derived from the values in the modular addition matrix,
  // where values contained in this matrix define the i-th t-base word in the integer to be reduced.
  // The following hard-coded method works for an arbitrary integer a and prime p = 2^64 - 2^8 - 1, where a < p^2.
  // Both a and p are represented as a collection of t-base words, where t is 2^8:
  // (a * b) = (ab_15, ab_14, ... , ab_1, ab_0) [16 bytes], p = (p_7, p_6, ... , p_1, p_0) [8 bytes]
  doubleword r = ((doubleword)(*a) * (doubleword)(*b));

  // Values (indices in value r, see modular addition matrix):
  // s0: 00 01 02 03 04 05 06 07
  // s1: 08 08 09 10 11 12 13 14
  // s2: 15 09 10 11 12 13 14 15
  // s3: 00 15 00 00 00 00 00 00

  word temp_var;
  doubleword x_tmp;
  x_tmp = r & 0xFFFFFFFFFFFFFFFF;
  /*
  //temp_var = ((r & 0x00000000000000FF0000000000000000) >> 64) | ((r & 0x00FFFFFFFFFFFFFF0000000000000000) >> 56);
  //temp_var = ((r & ((doubleword)(0xFF) << 64)) >> 64) | ((r & ((doubleword)(0x00FFFFFFFFFFFFFF) << 64)) >> 56);
  temp_var = ((r >> 64) & 0xFF) | (((r >> 64) & 0x00FFFFFFFFFFFFFF) << 8);
  x_tmp = x_tmp + temp_var;
  if(x_tmp >= BigIntLib::modulo_[0]) x_tmp -= BigIntLib::modulo_[0];
  //temp_var = ((r & 0xFF000000000000000000000000000000) >> 120) | ((r & 0xFFFFFFFFFFFFFF000000000000000000) >> 64);
  //temp_var = ((r & ((doubleword)(0xFF) << 120)) >> 120) | ((r & ((doubleword)(0x00FFFFFFFFFFFFFF) << 72)) >> 64);
  temp_var = ((r >> 120) & 0xFF) | (((r >> 72) & 0x00FFFFFFFFFFFFFF) << 8);
  x_tmp = x_tmp + temp_var;
  if(x_tmp >= BigIntLib::modulo_[0]) x_tmp -= BigIntLib::modulo_[0];
  //temp_var = ((r & 0xFF000000000000000000000000000000) >> 112);
  //temp_var = ((r & ((doubleword)(0xFF) << 120)) >> 112);
  temp_var = (((r >> 120) & 0xFF) << 8);
  x_tmp = x_tmp + temp_var;
  if(x_tmp >= BigIntLib::modulo_[0]) x_tmp -= BigIntLib::modulo_[0];
  *c = (word)x_tmp;
  */
  x_tmp = x_tmp + (((r >> 64) & 0xFF) | (((r >> 64) & 0x00FFFFFFFFFFFFFF) << 8)) +
    (((r >> 120) & 0xFF) | (((r >> 72) & 0x00FFFFFFFFFFFFFF) << 8)) +
    (((r >> 120) & 0xFF) << 8);
  if(x_tmp >= BigIntLib::modulo_[0]) x_tmp -= BigIntLib::modulo_[0];
  if(x_tmp >= BigIntLib::modulo_[0]) x_tmp -= BigIntLib::modulo_[0];
  if(x_tmp >= BigIntLib::modulo_[0]) x_tmp -= BigIntLib::modulo_[0];
  *c = (word)x_tmp;

  /*
  uchar* s = (uchar*)(&r);
  uint32 places = 8;
  uchar* s0;
  uchar s1[places];
  uchar* s2;
  uchar s3[8] = {0};
  
  s0 = s;
  memcpy(s1 + 1, s + 8, 7);
  s1[0] = s[8];
  s2 = s + 8;
  s2[0] = s[15];
  s3[1] = s[15];

  //ADD_FUNCTION(c, (word*)s0, (word*)s1);
  //ADD_FUNCTION(c, c, (word*)s2);
  //ADD_FUNCTION(c, c, (word*)s3);
  doubleword x = (doubleword)*((word*)s0) + *(word*)s1;
  if(x >= BigIntLib::modulo_[0]) x -= BigIntLib::modulo_[0];
  //if(x >= BigIntLib::modulo_[0]) x -= BigIntLib::modulo_[0];
  x = x + *(word*)s2;
  if(x >= BigIntLib::modulo_[0]) x -= BigIntLib::modulo_[0];
  //if(x >= BigIntLib::modulo_[0]) x -= BigIntLib::modulo_[0];
  x = x + *(word*)s3;
  if(x >= BigIntLib::modulo_[0]) x -= BigIntLib::modulo_[0];
  //if(x >= BigIntLib::modulo_[0]) x -= BigIntLib::modulo_[0];
  *c = (word)x;
  */
  
}

void BigIntLib::MulPF64FastCrandall(word* c, word* a, word* b) {
  doubleword r = ((doubleword)(*a) * (doubleword)(*b));
  word c_small = 59;
  doubleword q_i = r >> 64;
  r = r & 0xFFFFFFFFFFFFFFFF;
  doubleword t;

  // WHILE 1.1
  t = c_small * q_i;
  r = r + (t & 0xFFFFFFFFFFFFFFFF) * (q_i > 0);
  q_i = t >> 64;

  // WHILE 1.2
  t = c_small * q_i;
  r = r + (t & 0xFFFFFFFFFFFFFFFF) * (q_i > 0);

  // WHILE 2.1
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];

  // WHILE 2.2
  if(r >= BigIntLib::modulo_[0]) r -= BigIntLib::modulo_[0];

  *c = r;
}

void BigIntLib::AddPF136(word* c, word* a, word* b) {
  // Implement static ADD in PF with 136 bits...
  doubleword tmp;
  word carry = 0;
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
  //   tmp = (doubleword)a[i] + b[i] + carry;
  //   c[i] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  //   carry = (tmp >> WORD_SIZE);
  // }
  tmp = (doubleword)a[0] + b[0] + carry;
  c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);
  tmp = (doubleword)a[1] + b[1] + carry;
  c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);
  tmp = (doubleword)a[2] + b[2] + carry;
  c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);

  if(carry == 1 || BigIntLib::Greater(c, BigIntLib::modulo_)) {
    word borrow = 0;
    tmp = (doubleword)c[0] - BigIntLib::modulo_[0] - borrow;
    c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE) != 0;
    tmp = (doubleword)c[1] - BigIntLib::modulo_[1] - borrow;
    c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE) != 0;
    tmp = (doubleword)c[2] - BigIntLib::modulo_[2] - borrow;
    c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE) != 0;
  }
}

void BigIntLib::SubPF136(word* c, word* a, word* b) {
  // Implement static SUB in PF with 136 bits...
  doubleword tmp;
  word borrow = 0;
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
  //   tmp = (doubleword)a[i] - b[i] - borrow;
  //   c[i] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  //   borrow = (tmp >> WORD_SIZE) != 0;
  // }
  tmp = (doubleword)a[0] - b[0] - borrow;
  c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  tmp = (doubleword)a[1] - b[1] - borrow;
  c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  tmp = (doubleword)a[2] - b[2] - borrow;
  c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  
  if(borrow == 1) {
    word carry = 0;
    // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
    //   tmp = (doubleword)c[i] + BigIntLib::modulo_[i] + carry;
    //   c[i] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    //   carry = (tmp >> WORD_SIZE);
    // }
    tmp = (doubleword)c[0] + BigIntLib::modulo_[0] + carry;
    c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
    tmp = (doubleword)c[1] + BigIntLib::modulo_[1] + carry;
    c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
    tmp = (doubleword)c[2] + BigIntLib::modulo_[2] + carry;
    c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
  }
}

void BigIntLib::MulPF136Fast(word* c, word* a, word* b) {
  // Implement static MUL in PF with 136 bits followed by static reduction using a generalized mersenne prime...
  //memset(c, 0, BigIntLib::field_num_words_ * (WORD_SIZE / 8)); // Probably not needed
  word U;
  word V;
  doubleword UV;
  word c_temp[BigIntLib::field_num_words_ * 2];
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) c_temp[i] = 0;
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
  //   U = 0;
  //   for(uint32 j = 0; j < BigIntLib::field_num_words_; j++) {
  //     UV = c_temp[i + j] + ((doubleword)(a[i]) * (doubleword)(b[j])) + U;
  //     U = UV >> WORD_SIZE;
  //     V = (UV << WORD_SIZE) >> WORD_SIZE;
  //     c_temp[i + j] = V;
  //   }
  //   c_temp[i + BigIntLib::field_num_words_] = U;
  // }

  // STATIC
  memset(c_temp, 0, (BigIntLib::field_num_words_ * 2) * (WORD_SIZE / 8));
  // LOOP 1
  U = 0;
  UV = c_temp[0] + ((doubleword)(a[0]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[0] = V;
  UV = c_temp[1] + ((doubleword)(a[0]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[1] = V;
  UV = c_temp[2] + ((doubleword)(a[0]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[2] = V;
  c_temp[BigIntLib::field_num_words_] = U;
  // LOOP 2
  U = 0;
  UV = c_temp[1] + ((doubleword)(a[1]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[1] = V;
  UV = c_temp[2] + ((doubleword)(a[1]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[2] = V;
  UV = c_temp[3] + ((doubleword)(a[1]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[3] = V;
  c_temp[BigIntLib::field_num_words_ + 1] = U;
  // LOOP 3
  U = 0;
  UV = c_temp[2] + ((doubleword)(a[2]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[2] = V;
  UV = c_temp[3] + ((doubleword)(a[2]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[3] = V;
  UV = c_temp[4] + ((doubleword)(a[2]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[4] = V;
  c_temp[BigIntLib::field_num_words_ + 2] = U;

  // Modular additions
  // 17  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32
  // 33  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33
  // -1  33  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1

  uchar temp_var[24];
  uchar* pointer = (uchar*)c_temp;
  memset(temp_var, 0, 24);
  memset(c, 0, 24);
  memcpy(c, pointer, 17);
  memcpy(temp_var + 1, pointer + 17, 16);
  temp_var[0] = pointer[17];
  BigIntLib::Add(c, c, (word*)temp_var);
  memcpy(temp_var + 1, pointer + 18, 16);
  temp_var[0] = pointer[33];
  BigIntLib::Add(c, c, (word*)temp_var);
  memset(temp_var, 0, 24);
  temp_var[1] = pointer[33];
  BigIntLib::Add(c, c, (word*)temp_var);
}

void BigIntLib::AddPF256(word* c, word* a, word* b) {
  // Implement static ADD in PF with 256 bits...
  doubleword tmp;
  word carry = 0;
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
  //   tmp = (doubleword)a[i] + b[i] + carry;
  //   c[i] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  //   carry = (tmp >> WORD_SIZE);
  // }
  //BigIntLib::Print(a, 40);
  //BigIntLib::Print(b, 40);
  //BigIntLib::Print(c, 40);
  tmp = (doubleword)a[0] + b[0] + carry;
  c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);
  tmp = (doubleword)a[1] + b[1] + carry;
  c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);
  tmp = (doubleword)a[2] + b[2] + carry;
  c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);
  tmp = (doubleword)a[3] + b[3] + carry;
  c[3] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);

  if(carry == 1 || BigIntLib::Greater(c, BigIntLib::modulo_)) {
    word borrow = 0;
    // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
    //   tmp = (doubleword)c[i] - BigIntLib::modulo_[i] - borrow;
    //   c[i] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    //   borrow = (tmp >> WORD_SIZE) != 0;
    // }
    tmp = (doubleword)c[0] - BigIntLib::modulo_[0] - borrow;
    c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE) != 0;
    tmp = (doubleword)c[1] - BigIntLib::modulo_[1] - borrow;
    c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE) != 0;
    tmp = (doubleword)c[2] - BigIntLib::modulo_[2] - borrow;
    c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE) != 0;
    tmp = (doubleword)c[3] - BigIntLib::modulo_[3] - borrow;
    c[3] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE) != 0;
  }
}

void BigIntLib::SubPF256(word* c, word* a, word* b) {
  // Implement static SUB in PF with 256 bits...
  doubleword tmp;
  word borrow = 0;
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
  //   tmp = (doubleword)a[i] - b[i] - borrow;
  //   c[i] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  //   borrow = (tmp >> WORD_SIZE) != 0;
  // }
  tmp = (doubleword)a[0] - b[0] - borrow;
  c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  tmp = (doubleword)a[1] - b[1] - borrow;
  c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  tmp = (doubleword)a[2] - b[2] - borrow;
  c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  tmp = (doubleword)a[3] - b[3] - borrow;
  c[3] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  
  if(borrow == 1) {
    word carry = 0;
    // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
    //   tmp = (doubleword)c[i] + BigIntLib::modulo_[i] + carry;
    //   c[i] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    //   carry = (tmp >> WORD_SIZE);
    // }
    tmp = (doubleword)c[0] + BigIntLib::modulo_[0] + carry;
    c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
    tmp = (doubleword)c[1] + BigIntLib::modulo_[1] + carry;
    c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
    tmp = (doubleword)c[2] + BigIntLib::modulo_[2] + carry;
    c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
    tmp = (doubleword)c[3] + BigIntLib::modulo_[3] + carry;
    c[3] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
  }
}

void BigIntLib::MulPF256Fast(word* c, word* a, word* b) {
  // Implement static MUL in PF with 256 bits followed by static reduction using a generalized mersenne prime...
  //memset(c, 0, BigIntLib::field_num_words_ * (WORD_SIZE / 8)); // Probably not needed
  word U;
  word V;
  doubleword UV;
  word c_temp[BigIntLib::field_num_words_ * 2];
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) c_temp[i] = 0;
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
  //   U = 0;
  //   for(uint32 j = 0; j < BigIntLib::field_num_words_; j++) {
  //     UV = c_temp[i + j] + ((doubleword)(a[i]) * (doubleword)(b[j])) + U;
  //     U = UV >> WORD_SIZE;
  //     V = (UV << WORD_SIZE) >> WORD_SIZE;
  //     c_temp[i + j] = V;
  //   }
  //   c_temp[i + BigIntLib::field_num_words_] = U;
  // }

  // STATIC
  memset(c_temp, 0, (BigIntLib::field_num_words_ * 2) * (WORD_SIZE / 8));
  // LOOP 1
  U = 0;
  UV = c_temp[0] + ((doubleword)(a[0]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[0] = V;
  UV = c_temp[1] + ((doubleword)(a[0]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[1] = V;
  UV = c_temp[2] + ((doubleword)(a[0]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[2] = V;
  UV = c_temp[3] + ((doubleword)(a[0]) * (doubleword)(b[3])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[3] = V;
  c_temp[BigIntLib::field_num_words_] = U;
  // LOOP 2
  U = 0;
  UV = c_temp[1] + ((doubleword)(a[1]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[1] = V;
  UV = c_temp[2] + ((doubleword)(a[1]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[2] = V;
  UV = c_temp[3] + ((doubleword)(a[1]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[3] = V;
  UV = c_temp[4] + ((doubleword)(a[1]) * (doubleword)(b[3])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[4] = V;
  c_temp[BigIntLib::field_num_words_ + 1] = U;
  // LOOP 3
  U = 0;
  UV = c_temp[2] + ((doubleword)(a[2]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[2] = V;
  UV = c_temp[3] + ((doubleword)(a[2]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[3] = V;
  UV = c_temp[4] + ((doubleword)(a[2]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[4] = V;
  UV = c_temp[5] + ((doubleword)(a[2]) * (doubleword)(b[3])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[5] = V;
  c_temp[BigIntLib::field_num_words_ + 2] = U;
  // LOOP 4
  U = 0;
  UV = c_temp[3] + ((doubleword)(a[3]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[3] = V;
  UV = c_temp[4] + ((doubleword)(a[3]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[4] = V;
  UV = c_temp[5] + ((doubleword)(a[3]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[5] = V;
  UV = c_temp[6] + ((doubleword)(a[3]) * (doubleword)(b[3])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[6] = V;
  c_temp[BigIntLib::field_num_words_ + 3] = U;

  // Modular additions
  // 60  61  62  63  60  61  62  63  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  32  33  34  35  36  37  38  39  40
  // -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  41  42  43  44  45  46  47  48  49
  // -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  50  51  52  53  54  55  56  57  58
  // -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  59  60  61  62  63  -1  -1  -1  -1

  // Modular subtractions
  // 32  33  34  35  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59
  // 41  42  43  44  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63
  // 50  51  52  53  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  -1  -1  -1  -1  -1
  // 59  60  61  62  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  60  61  62  63  -1  -1  -1  -1  -1
  // -1  -1  -1  -1  50  51  52  53  54  55  56  57  58  59  60  61  62  63  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
  // -1  -1  -1  -1  54  55  56  57  58  59  60  61  62  63  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
  // -1  -1  -1  -1  59  60  61  62  63  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
  // -1  -1  -1  -1  63  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1

  // Following additions
  // This is slower than the optimized memcpy-versions found in the other field sizes,
  // but due to the comparatively high number of modular additions and subtractions it shouldn't make a huge difference.
  uchar temp_var[32];
  memcpy(c, c_temp, 32);
  uchar* pointer = (uchar*)c_temp;
  for(uint32 i = 0; i < BigIntLib::mod_addition_weight_; i++) {
    for(uint32 j = 0; j < 32; j++) { // 32 = 256 / 8, which is the word size for the reduction using the specified Generalized Mersenne Prime
      temp_var[j] = ((BigIntLib::mod_addition_matrix_[i][j]) == -1) ? 0 : pointer[(BigIntLib::mod_addition_matrix_[i][j])];
    }
    BigIntLib::Add(c, c, (word*)temp_var);
  }

  // Following subtractions
  for(uint32 i = 0; i < BigIntLib::mod_subtraction_weight_; i++) {
    for(uint32 j = 0; j < 32; j++) { // 32 = 256 / 8, which is the word size for the reduction using the specified Generalized Mersenne Prime
      temp_var[j] = ((BigIntLib::mod_subtraction_matrix_[i][j]) == -1) ? 0 : pointer[(BigIntLib::mod_subtraction_matrix_[i][j])];
    }
    BigIntLib::Sub(c, c, (word*)temp_var);
  }
  
  //memcpy(c, c_temp, 32);
}

void BigIntLib::AddPF272(word* c, word* a, word* b) {
  // Implement static ADD in PF with 272 bits...
  a[4] &= 0xFFFF; // Keep only first 16 bits of last word (16 + 4 * 64 = 272)
  b[4] &= 0xFFFF; // Keep only first 16 bits of last word (16 + 4 * 64 = 272)
  doubleword tmp;
  word carry = 0;
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
  //   tmp = (doubleword)a[i] + b[i] + carry;
  //   c[i] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  //   carry = (tmp >> WORD_SIZE);
  // }
  //BigIntLib::Print(a, 40);
  //BigIntLib::Print(b, 40);
  //BigIntLib::Print(c, 40);
  tmp = (doubleword)a[0] + b[0] + carry;
  c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);
  tmp = (doubleword)a[1] + b[1] + carry;
  c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);
  tmp = (doubleword)a[2] + b[2] + carry;
  c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);
  tmp = (doubleword)a[3] + b[3] + carry;
  c[3] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);
  tmp = (doubleword)a[4] + b[4] + carry;
  c[4] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  carry = (tmp >> WORD_SIZE);

  if(carry == 1 || BigIntLib::Greater(c, BigIntLib::modulo_)) { // Carry should never be 1 with an only partially used last word...
    word borrow = 0;
    // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
    //   tmp = (doubleword)c[i] - BigIntLib::modulo_[i] - borrow;
    //   c[i] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    //   borrow = (tmp >> WORD_SIZE) != 0;
    // }
    tmp = (doubleword)c[0] - BigIntLib::modulo_[0] - borrow;
    c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE) != 0;
    tmp = (doubleword)c[1] - BigIntLib::modulo_[1] - borrow;
    c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE) != 0;
    tmp = (doubleword)c[2] - BigIntLib::modulo_[2] - borrow;
    c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE) != 0;
    tmp = (doubleword)c[3] - BigIntLib::modulo_[3] - borrow;
    c[3] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE) != 0;
    tmp = (doubleword)c[4] - BigIntLib::modulo_[4] - borrow;
    c[4] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    borrow = (tmp >> WORD_SIZE);
  }
}

void BigIntLib::SubPF272(word* c, word* a, word* b) {
  // Implement static SUB in PF with 272 bits...
  a[4] &= 0xFFFF; // Keep only first 16 bits of last word (16 + 4 * 64 = 272)
  b[4] &= 0xFFFF; // Keep only first 16 bits of last word (16 + 4 * 64 = 272)
  doubleword tmp;
  word borrow = 0;
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
  //   tmp = (doubleword)a[i] - b[i] - borrow;
  //   c[i] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  //   borrow = (tmp >> WORD_SIZE) != 0;
  // }
  tmp = (doubleword)a[0] - b[0] - borrow;
  c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  tmp = (doubleword)a[1] - b[1] - borrow;
  c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  tmp = (doubleword)a[2] - b[2] - borrow;
  c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  tmp = (doubleword)a[3] - b[3] - borrow;
  c[3] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  tmp = (doubleword)a[4] - b[4] - borrow;
  c[4] = ((tmp << WORD_SIZE) >> WORD_SIZE);
  borrow = (tmp >> WORD_SIZE) != 0;
  
  if(borrow == 1) {
    word carry = 0;
    // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
    //   tmp = (doubleword)c[i] + BigIntLib::modulo_[i] + carry;
    //   c[i] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    //   carry = (tmp >> WORD_SIZE);
    // }
    tmp = (doubleword)c[0] + BigIntLib::modulo_[0] + carry;
    c[0] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
    tmp = (doubleword)c[1] + BigIntLib::modulo_[1] + carry;
    c[1] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
    tmp = (doubleword)c[2] + BigIntLib::modulo_[2] + carry;
    c[2] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
    tmp = (doubleword)c[3] + BigIntLib::modulo_[3] + carry;
    c[3] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
    tmp = (doubleword)c[4] + BigIntLib::modulo_[4] + carry;
    c[4] = ((tmp << WORD_SIZE) >> WORD_SIZE);
    carry = (tmp >> WORD_SIZE);
  }
}

void BigIntLib::MulPF272Fast(word* c, word* a, word* b) {
  // Implement static MUL in PF with 272 bits followed by static reduction using a generalized mersenne prime...
  //memset(c, 0, BigIntLib::field_num_words_ * (WORD_SIZE / 8)); // Probably not needed
  a[4] &= 0xFFFF; // Keep only first 16 bits of last word (16 + 4 * 64 = 272)
  b[4] &= 0xFFFF; // Keep only first 16 bits of last word (16 + 4 * 64 = 272)
  word U;
  word V;
  doubleword UV;
  word c_temp[BigIntLib::field_num_words_ * 2];
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) c_temp[i] = 0;
  // for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) {
  //   U = 0;
  //   for(uint32 j = 0; j < BigIntLib::field_num_words_; j++) {
  //     UV = c_temp[i + j] + ((doubleword)(a[i]) * (doubleword)(b[j])) + U;
  //     U = UV >> WORD_SIZE;
  //     V = (UV << WORD_SIZE) >> WORD_SIZE;
  //     c_temp[i + j] = V;
  //   }
  //   c_temp[i + BigIntLib::field_num_words_] = U;
  // }

  // STATIC
  memset(c_temp, 0, (BigIntLib::field_num_words_ * 2) * (WORD_SIZE / 8));
  // LOOP 1
  U = 0;
  UV = c_temp[0] + ((doubleword)(a[0]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[0] = V;
  UV = c_temp[1] + ((doubleword)(a[0]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[1] = V;
  UV = c_temp[2] + ((doubleword)(a[0]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[2] = V;
  UV = c_temp[3] + ((doubleword)(a[0]) * (doubleword)(b[3])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[3] = V;
  UV = c_temp[4] + ((doubleword)(a[0]) * (doubleword)(b[4])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[4] = V;
  c_temp[BigIntLib::field_num_words_] = U;
  // LOOP 2
  U = 0;
  UV = c_temp[1] + ((doubleword)(a[1]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[1] = V;
  UV = c_temp[2] + ((doubleword)(a[1]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[2] = V;
  UV = c_temp[3] + ((doubleword)(a[1]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[3] = V;
  UV = c_temp[4] + ((doubleword)(a[1]) * (doubleword)(b[3])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[4] = V;
  UV = c_temp[5] + ((doubleword)(a[1]) * (doubleword)(b[4])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[5] = V;
  c_temp[BigIntLib::field_num_words_ + 1] = U;
  // LOOP 3
  U = 0;
  UV = c_temp[2] + ((doubleword)(a[2]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[2] = V;
  UV = c_temp[3] + ((doubleword)(a[2]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[3] = V;
  UV = c_temp[4] + ((doubleword)(a[2]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[4] = V;
  UV = c_temp[5] + ((doubleword)(a[2]) * (doubleword)(b[3])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[5] = V;
  UV = c_temp[6] + ((doubleword)(a[2]) * (doubleword)(b[4])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[6] = V;
  c_temp[BigIntLib::field_num_words_ + 2] = U;
  // LOOP 4
  U = 0;
  UV = c_temp[3] + ((doubleword)(a[3]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[3] = V;
  UV = c_temp[4] + ((doubleword)(a[3]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[4] = V;
  UV = c_temp[5] + ((doubleword)(a[3]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[5] = V;
  UV = c_temp[6] + ((doubleword)(a[3]) * (doubleword)(b[3])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[6] = V;
  UV = c_temp[7] + ((doubleword)(a[3]) * (doubleword)(b[4])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[7] = V;
  c_temp[BigIntLib::field_num_words_ + 3] = U;
  // LOOP 5
  U = 0;
  UV = c_temp[4] + ((doubleword)(a[4]) * (doubleword)(b[0])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[4] = V;
  UV = c_temp[5] + ((doubleword)(a[4]) * (doubleword)(b[1])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[5] = V;
  UV = c_temp[6] + ((doubleword)(a[4]) * (doubleword)(b[2])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[6] = V;
  UV = c_temp[7] + ((doubleword)(a[4]) * (doubleword)(b[3])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[7] = V;
  UV = c_temp[8] + ((doubleword)(a[4]) * (doubleword)(b[4])) + U;
  U = UV >> WORD_SIZE;
  V = (UV << WORD_SIZE) >> WORD_SIZE;
  c_temp[8] = V;
  c_temp[BigIntLib::field_num_words_ + 4] = U;

  // Reduction
  uchar* pointer = (uchar*)c_temp;
  word temp_var[5] = {0}; // 34 = 272 / 8, which is the word size for the reduction using the specified Generalized Mersenne Prime
  uchar* temp_var_pointer = (uchar*)temp_var;
  for(uint32 i = 0; i < BigIntLib::field_num_words_; i++) { // better with memset
    c[i] = 0;
    temp_var[i] = 0;
  }

  // Modular Additions
  // 34  35  36  37  38  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62
  // 63  64  65  66  67  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67
  // -1  -1  -1  -1  -1  63  64  65  66  67  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1
  
  // c <- s0
  memcpy(c, pointer, 34);

  // First addition
  memcpy(temp_var_pointer, pointer + 34, 5);
  memcpy(temp_var_pointer + 5, pointer + 34, 29);
  BigIntLib::Add(c, c, temp_var);

  // Second addition
  memcpy(temp_var_pointer, pointer + 63, 5);
  memcpy(temp_var_pointer + 5, pointer + 39, 29);
  BigIntLib::Add(c, c, temp_var);

  // Third addition
  memset(temp_var_pointer, 0, 34);
  memcpy(temp_var_pointer + 5, pointer + 63, 5);
  BigIntLib::Add(c, c, temp_var);

  // Following additions
  // uchar temp_v[34];
  // memcpy(c, c_temp, BigIntLib::field_size_bytes_);
  // uchar* ptr = (uchar*)c_temp;
  // for(uint32 i = 0; i < BigIntLib::mod_addition_weight_; i++) {
  //   for(uint32 j = 0; j < 34; j++) { // 34 = 272 / 8, which is the word size for the reduction using the specified Generalized Mersenne Prime
  //     temp_v[j] = ((BigIntLib::mod_addition_matrix_[i][j]) == -1) ? 0 : ptr[(BigIntLib::mod_addition_matrix_[i][j])];
  //   }
  //   BigIntLib::Add(c, c, (word*)temp_v);
  // }

  //memcpy(c, c_temp, BigIntLib::field_size_bytes_);
}

void BigIntLib::Times2PF(word* c, word* a) {
  BigIntLib::Add(c, a, a);
}

void BigIntLib::Times3PF(word* c, word* a) {
  BigIntLib::Add(c, a, a);
  BigIntLib::Add(c, c, a);
}

void BigIntLib::XorBF3(word* c, word* a, word* b) {
  *c = *a ^ *b;
}

void BigIntLib::MulBF3Fast(word* c, word* a, word* b) {
  /*
  doubleword r = b[0];
  asm("pclmulqdq %2, %1, %0;"
    : "+x"(r)
    : "x"(a[0]), "i"(0)
    );
  */
  doubleword r = 0;
  asm("pclmulqdq %2, %1, %0;"
    : "=x"(r)
    : "x"(a[0]), "i"(0), "0"(b[0])
    );
  
  word c0 = r & 0x7; // LS 3 bits
  word c1 = (r & 0x18) >> 3; // MS 2 bits

  // Add c1 to bits 1 and 0 of c0 (c1 doesn't neet to be done first, because it's not affected)
  *c = c0 ^ (c1 << 1) ^ c1;
}

void BigIntLib::Times2BF3(word* c, word* a) {
  word r = 0;
  memcpy(&r, a, BigIntLib::field_size_bytes_);
  r <<= 1;
  // Bits 1, 0
  word T = r >> 3;
  *c = (r ^ (T << 1) ^ T) & 0x7;
}

void BigIntLib::Times3BF3(word* c, word* a) {
  word r = 0;
  memcpy(&r, a, BigIntLib::field_size_bytes_);
  r <<= 2;
  // Bits 1, 0
  word T = r >> 3;
  *c = (r ^ (T << 1) ^ T) & 0x7;
}

void BigIntLib::XorBF17(word* c, word* a, word* b) {
  *c = *a ^ *b;
}

void BigIntLib::MulBF17Fast(word* c, word* a, word* b) {
  
  doubleword r = 0;
  asm("pclmulqdq %2, %1, %0;"
    : "=x"(r)
    : "x"(a[0]), "i"(0), "0"(b[0])
    );
  
  // Reduction
  // p(x) = x^17 + x^3 + 1
  /*
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
  *c = r & 0x1FFFF;
  */
  word c0 = r & 0x1FFFF; // LS 17 bits
  word c1 = r >> 17; // MS 16 bits
  
  c1 = c1 ^ (c1 >> 14);
  // Add c1 to bits 3 and 0 of c0
  c0 = c0 ^ (c1 << 3) ^ c1;
  
  // Build result
  *c = c0 & 0x1FFFF;
}

void BigIntLib::Times2BF17(word* c, word* a) {
  word r = 0;
  memcpy(&r, a, BigIntLib::field_size_bytes_);
  r <<= 1;
  // Bits 3, 0
  word T = r >> 17;
  *c = (r ^ (T << 3) ^ T) & 0x1FFFF;
}

void BigIntLib::Times3BF17(word* c, word* a) {
  word r = 0;
  memcpy(&r, a, BigIntLib::field_size_bytes_);
  r <<= 2;
  // Bits 3, 0
  word T = r >> 17;
  *c = (r ^ (T << 3) ^ T) & 0x1FFFF;
}

void BigIntLib::XorBF33(word* c, word* a, word* b) {
  *c = *a ^ *b;
}

void BigIntLib::MulBF33Fast(word* c, word* a, word* b) {
  /*
  doubleword r = 0;
  word b_temp = *b;

  // Fast multiplication (Left-to-right comb method with windows of width w)
  // Precompute u(x) * b(x) for all possible values of u(x) (w = 4, 16 possible values)
  // This takes (2 * 16) - 15 = 17 XOR operations and (3 * 16) / 2 = 24 shifts
  // Note: This can be implemented faster by reusing already computed results! E.g. B_u[3] = b_temp ^ (b_temp << 1) = b_temp ^ B_u[2]
  doubleword B_u[16];
  // 0: 0b0000
  B_u[0] = 0x0;
  // 1: 0b0001
  B_u[1] = b_temp;
  // 2: 0b0010
  B_u[2] = (b_temp << 1);
  // 3: 0b0011
  B_u[3] = b_temp ^ (b_temp << 1);
  // 4: 0b0100
  B_u[4] = (b_temp << 2);
  // 5: 0b0101
  B_u[5] = b_temp ^ (b_temp << 2);
  // 6: 0b0110
  B_u[6] = (b_temp << 1) ^ (b_temp << 2);
  // 7: 0b0111
  B_u[7] = b_temp ^ (b_temp << 1) ^ (b_temp << 2);
  // 8: 0b1000
  B_u[8] = (b_temp << 3);
  // 9: 0b1001
  B_u[9] = b_temp ^ (b_temp << 3);
  // 10: 0b1010
  B_u[10] = (b_temp << 1) ^ (b_temp << 3);
  // 11: 0b1011
  B_u[11] = b_temp ^ (b_temp << 1) ^ (b_temp << 3);
  // 12: 0b1100
  B_u[12] = (b_temp << 2) ^ (b_temp << 3);
  // 13: 0b1101
  B_u[13] = b_temp ^ (b_temp << 2) ^ (b_temp << 3);
  // 14: 0b1110
  B_u[14] = (b_temp << 1) ^ (b_temp << 2) ^ (b_temp << 3);
  // 15: 0b1111
  B_u[15] = b_temp ^ (b_temp << 1) ^ (b_temp << 2) ^ (b_temp << 3);
  */

  // Unrolled version of the multiplication loop
  /*
  r = r ^ B_u[(*a >> 32) & 0xF];
  r = r << 4;
  r = r ^ B_u[(*a >> 28) & 0xF];
  r = r << 4;
  r = r ^ B_u[(*a >> 24) & 0xF];
  r = r << 4;
  r = r ^ B_u[(*a >> 20) & 0xF];
  r = r << 4;
  r = r ^ B_u[(*a >> 16) & 0xF];
  r = r << 4;
  r = r ^ B_u[(*a >> 12) & 0xF];
  r = r << 4;
  r = r ^ B_u[(*a >> 8) & 0xF];
  r = r << 4;
  r = r ^ B_u[(*a >> 4) & 0xF];
  r = r << 4;
  r = r ^ B_u[(*a) & 0xF];
  */
  
  /*
  // Optimized version
  r = ((((((((((((((((B_u[(*a >> 32) & 0xF]) << 4) ^
    B_u[(*a >> 28) & 0xF]) << 4) ^
    B_u[(*a >> 24) & 0xF]) << 4) ^
    B_u[(*a >> 20) & 0xF]) << 4) ^
    B_u[(*a >> 16) & 0xF]) << 4) ^
    B_u[(*a >> 12) & 0xF]) << 4) ^
    B_u[(*a >> 8) & 0xF]) << 4) ^
    B_u[(*a >> 4) & 0xF]) << 4) ^
    B_u[(*a) & 0xF];
  */

  // (Experimental) Carry-Less multiplication using CPU instruction PCLMULQDQ
  /*
  doubleword r = b[0];
  asm("pclmulqdq %2, %1, %0;"
    : "+x"(r)
    : "x"(a[0]), "i"(0)
    );
  */
  doubleword r = 0;
  asm("pclmulqdq %2, %1, %0;"
    : "=x"(r)
    : "x"(a[0]), "i"(0), "0"(b[0])
    );

  
  
  // Reduction
  // x^33 = x^6 + x^3 + x + 1 mod p(x)
  // x^64 = x^37 + x^34 + x^32 + x^31 mod p(x)
  // Manipulate from higher words to lower words
  
  word c0 = r & 0x1FFFFFFFF; // LS 33 bits
  word c1 = r >> 33; // MS 32 bits

  word T = c1;
  c1 = c1 ^ (T >> 27) ^ (T >> 30); // 27 = 33 - 6, 30 = 33 - 3, 32 = 33 - 1 (omitted, all zeros), x^0 does not affect c1
  T = c1;
  *c = c0 ^ ((T << 6) & 0x1FFFFFFFF) ^ ((T << 3) & 0x1FFFFFFFF) ^ ((T << 1) & 0x1FFFFFFFF) ^ T; // = c0, for x^6, x^3, x^1, x^0

  // p(x) = x^33 + x^10 + 1
  /*
  word r_0_n_1 = r & 0x1FFFFFFFF;
  word r_n_deg_r = r >> 33;
  word x_m_1 = 0x401; // = 2^10 + 1
  doubleword t = 0;
  asm("pclmulqdq %2, %1, %0;"
    : "=x"(t)
    : "x"(r_n_deg_r), "i"(0), "0"(x_m_1)
    );
  r = r_0_n_1 ^ t;
  t = r >> 33;
  r = r ^ (t << 10) ^ t;
  *c = r & 0x1FFFFFFFF;
  */

  /*
  uint32 c0 = r & 0xFFFFFFFF; // first 32 bits
  uint32 c1 = (r >> 32) & 0xFFFFFFFF; // second 32 bits
  uint32 c2 = (r >> 64) & 0xFFFFFFFF; // third 32 bits (only 1 bit)

  // Reduce LSB of c2 and add to bits 37, 34, 32 and 31 of c (because x^64 = x^37 + x^34 + x^32 + x^31 mod p(x))
  uint32 T = c2;
  c1 = c1 ^ (T << 5);
  c1 = c1 ^ (T << 2);
  c1 = c1 ^ T;
  c0 = c0 ^ (T << 31);

  // Reduce MS 31 bits of c1 and add to bits 6, 3, 1 and 0 of c (because x^33 = x^6 + x^3 + x + 1 mod p(x))
  // Parts relevant for c1
  T = c1 >> 1;
  c1 = c1 ^ (T >> 26);
  c1 = c1 ^ (T >> 29);
  c1 = c1 ^ (T >> 31);
  // Parts relevant for c0
  T = c1 >> 1;
  c0 = c0 ^ (T << 6);
  c0 = c0 ^ (T << 3);
  c0 = c0 ^ (T << 1);
  c0 = c0 ^ T;

  // Build result = (LSB of c1 || c0)
  *c = ((word)(c1 & 0x1) << 32) | c0;
  */
}

void BigIntLib::Times2BF33(word* c, word* a) {
  word r = 0;
  memcpy(&r, a, BigIntLib::field_size_bytes_);
  r <<= 1;
  // Bits 6, 3, 1, 0
  word T = r >> 33;
  *c = (r ^ (T << 6) ^ (T << 3) ^ (T << 1) ^ T) & 0x1FFFFFFFF;
}

void BigIntLib::Times3BF33(word* c, word* a) {
  word r = 0;
  memcpy(&r, a, BigIntLib::field_size_bytes_);
  r <<= 2;
  // Bits 6, 3, 1, 0
  word T = r >> 33;
  *c = (r ^ (T << 6) ^ (T << 3) ^ (T << 1) ^ T) & 0x1FFFFFFFF;
}

void BigIntLib::XorBF65(word* c, word* a, word* b) {
  c[0] = a[0] ^ b[0];
  c[1] = a[1] ^ b[1];
}

void BigIntLib::MulBF65Fast(word* c, word* a, word* b) {
  //word r[3] = {0};
  //doubleword b_temp = *((doubleword*)b);

  // Slow multiplication (Right-to-left comb method)
  /*
  for(uint32 i = 0; i < 64; i++) { // Make this operation faster with window-based multiplication
    if((((*a) & ((word)0x1 << i)) >> i) == 0x1)
      r = r ^ b_temp;
    if(i != 63)
      b_temp = b_temp << 1;
  }
  */
  

  // Fast multiplication (Left-to-right comb method with windows of width w)
  // Precompute u(x) * b(x) for all possible values of u(x) (w = 4, 16 possible values)
  // This takes (2 * 16) - 15 = 17 XOR operations and (3 * 16) / 2 = 24 shifts
  // Note: This can be implemented faster by reusing already computed results! E.g. B_u[3] = b_temp ^ (b_temp << 1) = b_temp ^ B_u[2]
  /*
  doubleword B_u[16];
  // 0: 0b0000
  B_u[0] = 0x0;
  // 1: 0b0001
  B_u[1] = b_temp;
  // 2: 0b0010
  B_u[2] = (b_temp << 1);
  // 3: 0b0011
  B_u[3] = b_temp ^ (b_temp << 1);
  // 4: 0b0100
  B_u[4] = (b_temp << 2);
  // 5: 0b0101
  B_u[5] = b_temp ^ (b_temp << 2);
  // 6: 0b0110
  B_u[6] = (b_temp << 1) ^ (b_temp << 2);
  // 7: 0b0111
  B_u[7] = b_temp ^ (b_temp << 1) ^ (b_temp << 2);
  // 8: 0b1000
  B_u[8] = (b_temp << 3);
  // 9: 0b1001
  B_u[9] = b_temp ^ (b_temp << 3);
  // 10: 0b1010
  B_u[10] = (b_temp << 1) ^ (b_temp << 3);
  // 11: 0b1011
  B_u[11] = b_temp ^ (b_temp << 1) ^ (b_temp << 3);
  // 12: 0b1100
  B_u[12] = (b_temp << 2) ^ (b_temp << 3);
  // 13: 0b1101
  B_u[13] = b_temp ^ (b_temp << 2) ^ (b_temp << 3);
  // 14: 0b1110
  B_u[14] = (b_temp << 1) ^ (b_temp << 2) ^ (b_temp << 3);
  // 15: 0b1111
  B_u[15] = b_temp ^ (b_temp << 1) ^ (b_temp << 2) ^ (b_temp << 3);
  */

  // Compute for each w-bit window, where w = 4
  /*
  uint32 u = 0;
  for(int i = 15; i >= 0; i--) { // 16 = W / w = 64 / 4
    u = (*a >> (i * 4)) & 0xF;
    r = r ^ B_u[u];
    if(i != 0)
      r = r << 4;
  }
  */
  // Unrolled version of the multiplication loop
  //doubleword* r_dw_pointer = (doubleword*)r;
  /*
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[1]) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 60) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 56) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 52) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 48) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 44) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 40) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 36) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 32) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 28) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 24) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 20) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 16) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 12) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 8) & 0xF];
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0] >> 4) & 0xF];
  */

  // Optimized version
  /*
  doubleword* r_dw_pointer = (doubleword*)r;
  *r_dw_pointer = ((((((((((((((((((((((((((((((B_u[(a[1]) & 0xF]) << 4) ^
    B_u[(a[0] >> 60) & 0xF]) << 4) ^
    B_u[(a[0] >> 56) & 0xF]) << 4) ^
    B_u[(a[0] >> 52) & 0xF]) << 4) ^
    B_u[(a[0] >> 48) & 0xF]) << 4) ^
    B_u[(a[0] >> 44) & 0xF]) << 4) ^
    B_u[(a[0] >> 40) & 0xF]) << 4) ^
    B_u[(a[0] >> 36) & 0xF]) << 4) ^
    B_u[(a[0] >> 32) & 0xF]) << 4) ^
    B_u[(a[0] >> 28) & 0xF]) << 4) ^
    B_u[(a[0] >> 24) & 0xF]) << 4) ^
    B_u[(a[0] >> 20) & 0xF]) << 4) ^
    B_u[(a[0] >> 16) & 0xF]) << 4) ^
    B_u[(a[0] >> 12) & 0xF]) << 4) ^
    B_u[(a[0] >> 8) & 0xF]) << 4) ^
    B_u[(a[0] >> 4) & 0xF];

  // Modified << 4 (store MS 4 bits of current value to r[2])
  r[2] = (*r_dw_pointer >> 124) & 0xF;
  *r_dw_pointer = *r_dw_pointer << 4;
  *r_dw_pointer = *r_dw_pointer ^ B_u[(a[0]) & 0xF];
  */

  // Word-wise polynomial multiplication using the CLMUL instruction
  word r[3];
  //memset(r, 0x0, 3 * 8);
  /*
  doubleword temp_res = a[0];
  asm("pclmulqdq %2, %1, %0;"
    : "+x"(temp_res)
    : "x"(b[0]), "i"(0)
    );
  */

  doubleword temp_res = 0;
  asm("pclmulqdq %2, %1, %0;"
    : "=x"(temp_res)
    : "x"(a[0]), "i"(0), "0"(b[0])
    );


  r[0] = ((temp_res << 64) >> 64);
  r[1] = (temp_res >> 64);
  //temp_res = a[0] * b[1];
  //r[1] = r[1] ^ ((temp_res << 64) >> 64); // or: r[1] = r[1] ^ (a[0] * b[1]), since b[1] is either 0 or 1
  r[1] = r[1] ^ (a[0] * b[1]);
  //temp_res = a[1] * b[0];
  //r[1] = r[1] ^ ((temp_res << 64) >> 64); // or: r[1] = r[1] ^ (a[1] * b[0]), since a[1] is either 0 or 1
  r[1] = r[1] ^ (a[1] * b[0]);
  //temp_res = a[1] * b[1];
  //r[2] = ((temp_res << 64) >> 64); // or: r[2] = r[2] ^ (a[1] * b[1]), since both a[1] and b[1] are either 0 or 1
  r[2] = a[1] * b[1];
  //std::cout << "a = 0x" << BigIntLib::ToString(a, 9) << std::endl;
  //std::cout << "b = 0x" << BigIntLib::ToString(b, 9) << std::endl;
  //std::cout << "r = " << BigIntLib::ToString(r, 17) << std::endl;

  // Reduction
  // x^65 = x^4 + x^3 + x + 1 mod p(x)
  // x^128 = x^67 + x^66 + x^64 + x^63 mod p(x)
  // Manipulate from higher words to lower words

  //std::cout << "r: " << BigIntLib::ToString(&r, BigIntLib::field_size_bytes_ * 2) << std::endl;
  
  /*
  doubleword c0 = ((doubleword)(r[1] & 0x1) << 64) | r[0]; // LS 65 bits
  word c1 = ((r[2] & 0x1) << 63) | (r[1] >> 1); // MS 64 bits

  doubleword T = c1;
  c1 = c1 ^ (T >> 61) ^ (T >> 62); // 61 = 65 - 4, 62 = 65 - 3, 64 = 65 - 1 (omitted, all zeros), x^0 does not affect c1
  T = c1;
  c0 = c0 ^ ((T << 63) >> 59) ^ ((T << 63) >> 60) ^ ((T << 63) >> 62) ^ T; // for x^4, x^3, x^1, x^0
  c[0] = c0 & 0xFFFFFFFFFFFFFFFF; // LS 64 bits of c0
  c[1] = (c0 >> 64) & 0x1; // 64th bit of c0 (x^64)
  */
  
  word c0 = r[0]; // first 64 bits
  word c1 = r[1]; // second 64 bits
  word c2 = r[2]; // third 64 bits (only 1 bit)

  // Reduce LSB of c2 and add to bits 67, 66, 64 and 63 of c (because x^128 = x^67 + x^66 + x^64 + x^63 mod p(x))
  word T = c2;
  c1 = c1 ^ (T << 3);
  c1 = c1 ^ (T << 2);
  c1 = c1 ^ T;
  c0 = c0 ^ (T << 63);

  // Reduce MS 63 bits of c1 and add to bits 4, 3, 1 and 0 of c (because x^65 = x^4 + x^3 + x + 1 mod p(x))
  // Parts relevant for c0
  T = c1 >> 1;
  c1 = c1 ^ (T >> 60);
  c1 = c1 ^ (T >> 61);
  c1 = c1 ^ (T >> 63);
  // Parts relevant for c0
  T = c1 >> 1;
  c0 = c0 ^ (T << 4);
  c0 = c0 ^ (T << 3);
  c0 = c0 ^ (T << 1);
  c0 = c0 ^ T;

  // Build result = (LSB of c1 || c0)
  c[1] = c1 & 0x1;
  c[0] = c0;

  /*
  // Slower anyway, small bug inside
  doubleword red_bits = (r[1] >> 1) | (r[2] << 63);
  doubleword* r_val = (doubleword*)r;
  *r_val = *r_val ^ (red_bits << 4);
  red_bits = (r[1] >> 1) | (r[2] << 63);
  *r_val = *r_val ^ (red_bits << 3);
  red_bits = (r[1] >> 1) | (r[2] << 63);
  *r_val = *r_val ^ (red_bits << 1);
  red_bits = (r[1] >> 1) | (r[2] << 63);
  *r_val = *r_val ^ (red_bits);

  c[1] = r[1] & 0x1;
  c[0] = r[0];
  */
  
}

void BigIntLib::Times2BF65(word* c, word* a) {
  doubleword r = 0;
  memcpy(&r, a, BigIntLib::field_size_bytes_);
  r <<= 1;
  // Bits 4, 3, 1, 0
  doubleword T = r >> 65;
  *c = ((r ^ (T << 4) ^ (T << 3) ^ (T << 1) ^ T) << 63) >> 63;
}

void BigIntLib::Times3BF65(word* c, word* a) {
  doubleword r = 0;
  memcpy(&r, a, BigIntLib::field_size_bytes_);
  r <<= 2;
  // Bits 4, 3, 1, 0
  doubleword T = r >> 65;
  *c = ((r ^ (T << 4) ^ (T << 3) ^ (T << 1) ^ T) << 63) >> 63;
}
