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

// Functions for operations in prime fields and binary fields
#ifndef BIGINTLIB_H
#define BIGINTLIB_H

#include <chrono>
#include <iostream>
#include <openssl/evp.h>

#define WORD_SIZE 64

typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned long uint64;
typedef __uint128_t uint128;
typedef unsigned char uchar;
typedef uint64 word;
typedef uint128 doubleword;

class BigIntLib {
public:
  static void Init(uint32 branch_size_bits, uint32 field_type);
  static void CleanUp();

  // Wrappers
  static void Add(word* c, word* a, word* b);
  static void AddSpec(word* c, word* a, word* b); // Addition without reduction
  static void Sub(word* c, word* a, word* b);
  static void Mul(word* c, word* a, word* b);
  static void MulSpec(word* c, word* a); // Multiplication with a constant value
  static void SolinasReduc(word* out, word* in);
  static void DoubleAdd(word* c, word* a, word* b);
  static void Times2(word* c, word* a);
  static void Times3(word* c, word* a);
  static bool IsZero(word* a);
  static bool Smaller(word* a, word* b); // Not yet a wrapper, probably not needed as wrapper, neither is Greater(.)
  static bool Greater(word* a, word* b);
  static uchar Compare(word* a, word* b); // Not yet a wrapper, probably not needed as wrapper, neither is Greater(.)
  static void TryReduce(word* a);
  static void GetRandomFieldElements(void* destination, uchar* seed, uint32 num_elements);
  static void GetRandomFieldElementsNonZero(void* destination, uchar* seed, uint32 num_elements);

  // Util
  static void SetSeed(void* source, uint32 num_bytes);
  static void FillRandom(void* destination, uchar* seed, uint32 num_bytes);
  static void Print(void* source, uint32 num_bytes);
  static std::string ToString(void* source, uint32 num_bytes = 0); // Always in field size

  // Variables
  static uint32 field_type_;
  static uint32 field_size_bits_;
  static uint32 field_size_bytes_;
  static uint32 field_num_words_;
  static uint64 msb_word_mask_;
  static uint32 fastreduc_num_words_;
  static word* modulo_;
  static EVP_CIPHER_CTX* evp_cipher_ctx_;
  static uchar* evp_cipher_buffer_;
  static int** reduction_matrix_;
  static int** reduction_matrix_positive_;
  static int** reduction_matrix_negative_;
  static uint32 mod_addition_weight_;
  static int** mod_addition_matrix_;
  static uint32 mod_subtraction_weight_;
  static int** mod_subtraction_matrix_;

private:
  BigIntLib();
  ~BigIntLib();

  // Helpers for reduction
  static void ComputeReductionMatrix(int* t_coeff);
  static void ComputeModularAdditionMatrix();
  static void ComputeModularSubtractionMatrix();

  // Utils
  static bool GreaterPF(word* a, word* b);
  static bool GreaterBF(word* a, word* b);
  static void TryReducePF(word* a);
  static void TryReduceBF(word* a);
  static void GetRandomFieldElementsPF(void* destination, uchar* seed, uint32 num_elements);
  static void GetRandomFieldElementsBF(void* destination, uchar* seed, uint32 num_elements);
  static void GetRandomFieldElementsNonZeroPF(void* destination, uchar* seed, uint32 num_elements); // Probably slightly redundant
  static void GetRandomFieldElementsNonZeroBF(void* destination, uchar* seed, uint32 num_elements); // Probably slightly redundant

  // Prime field with 3 bits
  static void AddPF3(word* c, word* a, word* b);
  static void SubPF3(word* c, word* a, word* b);
  static void MulPF3Fast(word* c, word* a, word* b);

  // Prime field with 4 bits
  static void AddPF4(word* c, word* a, word* b);
  static void SubPF4(word* c, word* a, word* b);
  static void MulPF4Fast(word* c, word* a, word* b);
  static void MulPF4FastCrandall(word* c, word* a, word* b);
  static void DoubleAddPF4(word* c, word* a, word* b);

  // Prime field with 16 bits
  static void AddPF16(word* c, word* a, word* b);
  static void AddSpecPF16(word* c, word* a, word* b);
  static void SubPF16(word* c, word* a, word* b);
  static void MulPF16Fast(word* c, word* a, word* b);
  static void MulSpecPF16(word* c, word* a);
  static void SolinasReducPF16(word* out, word* in);
  static void MulPF16FastCrandall(word* c, word* a, word* b);

  // Prime field with 32 bits
  static void AddPF32(word* c, word* a, word* b);
  static void SubPF32(word* c, word* a, word* b);
  static void MulPF32Fast(word* c, word* a, word* b);
  static void MulPF32FastCrandall(word* c, word* a, word* b);
  
  // Prime field with 64 bits
  static void AddPF64(word* c, word* a, word* b);
  static void SubPF64(word* c, word* a, word* b);
  static void MulPF64(word* c, word* a, word* b);
  static void MulPF64Fast(word* c, word* a, word* b);
  static void MulPF64FastCrandall(word* c, word* a, word* b);

  // Prime field with 136 bits
  static void AddPF136(word* c, word* a, word* b);
  static void SubPF136(word* c, word* a, word* b);
  static void MulPF136Fast(word* c, word* a, word* b);

  // Prime field with 256 bits
  static void AddPF256(word* c, word* a, word* b);
  static void SubPF256(word* c, word* a, word* b);
  static void MulPF256Fast(word* c, word* a, word* b);

  // Prime field with 272 bits
  static void AddPF272(word* c, word* a, word* b);
  static void SubPF272(word* c, word* a, word* b);
  static void MulPF272Fast(word* c, word* a, word* b);

  // Prime fields
  static void Times2PF(word* c, word* a);
  static void Times3PF(word* c, word* a);

  // Binary field with 3 bits
  static void XorBF3(word* c, word* a, word* b);
  static void MulBF3Fast(word* c, word* a, word* b);
  static void Times2BF3(word* c, word* a);
  static void Times3BF3(word* c, word* a);

  // Binary field with 17 bits
  static void XorBF17(word* c, word* a, word* b);
  static void MulBF17Fast(word* c, word* a, word* b);
  static void Times2BF17(word* c, word* a);
  static void Times3BF17(word* c, word* a);

  // Binary field with 33 bits
  static void XorBF33(word* c, word* a, word* b);
  static void MulBF33Fast(word* c, word* a, word* b);
  static void Times2BF33(word* c, word* a);
  static void Times3BF33(word* c, word* a);

  // Binary field with 65 bits
  static void XorBF65(word* c, word* a, word* b);
  static void MulBF65Fast(word* c, word* a, word* b);
  static void Times2BF65(word* c, word* a);
  static void Times3BF65(word* c, word* a);
};

#endif // BIGINTLIB_H