// Common includes and definitions between files
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

#ifndef COMMON_H
#define COMMON_H

#include <chrono>
#include <iostream>

typedef unsigned int uint32;
typedef unsigned long uint64;
typedef unsigned char uchar;

// This struct contains all values of one view (branch)
typedef struct {
  uchar* values_;
  uint32 num_gates_;
  uint32 gate_size_;
} View;

// This struct is returned from the Circuit to ZKBPP::sign(.)
typedef struct {
  uchar* x_3_;
  uchar** key_shares_;
  uchar* random_tapes_;
  uchar* random_tapes_hashs_;
  uchar** views_;
  //uchar* views_;
  uchar* y_;
  uchar* y_shares_;
} SignData;

// This struct is returned from the Circuit to ZKBPP:verify(.)
typedef struct {
  uchar* view_;
  uchar* y_share_; // y share of the given branch
  uchar* y_e2_; // Has to be calculated from circuit, y_e+2 = y + y_e + y_e+1
} VerifyData;

// This struct contains all SignData* from the iterations
typedef struct {
  SignData** sds_;
} ContainerSignData;

// This struct contains all VerifyData* from the iterations
typedef struct {
  VerifyData** vds_;
} ContainerVerifyData;

// C_ij = [H'(k_ij, View_ij)]
typedef struct {
  uchar* H_k_View_;
} C;

// D_ij = [k_ij || View_ij]
typedef struct {
  uchar* k_View_;
} D;

// Container for C's and D's (all iterations, C_ij and D_ij)
typedef struct {
  C*** Cs_; // Matrix, because (party_size * iterations) entries
  D*** Ds_; // Matrix, because (party_size * iterations) entries
} ContainerCD;

// a_ij = [y_i1, y_i2, y_i3, C_i1, C_i2, C_i3] (concatenation)
typedef struct {
  uchar* ys_C_hashs_;
} A;

// Container for A structs
typedef struct {
  A** as_;
} ContainerA;

// b_i = [y_(e_i+2), C_(e_i+2)], where the value from C ist stored directly
typedef struct {
  uchar* y_e2_;
  uchar* H_k_View_;
} B;

// z_i = [View_i+1, k_i, k_i+1, x_3] where e = i, + key_shares + y share
typedef struct {
  uchar* view_;
  uchar* k_1_;
  uchar* k_2_;
  uchar* k_1_hash_;
  uchar* k_2_hash_;
  uchar* x_3_; // Not used if e_i = 0
  uchar** key_shares_;
  uchar* y_share_;
} Z;

// This struct is returned from ZKBOO to the user (and can be used with ZKBPP::Verify(.))
// Proof contains the challenge E, all b_i's and all z_i's (i for iteration)
// p = [e, (b_1, z_1), ..., (b_t, z_t)]
typedef struct {
  uint32 num_iterations_;
  uchar* e_;
  // b_i Part
  uchar** y_e2_;
  uchar** H_k_View_;
  //B** bs_;
  Z** zs_;
  uchar** ys_; // Added, not necessary, only for testing (should all be the same)
} Proof;

// Forward declarations
class Circuit;

#endif // COMMON_H