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

#include "ZKBPP.h"
#include "BigIntLib.h"
#include <cstring>
#include <cmath> // log2(.)
#include <iomanip> // For HEX output
#include <thread>
#ifndef __APPLE__
#include <omp.h> // Iteration parallelization
#endif
#include "openssl/sha.h" // For SHA-256
#include <openssl/rand.h> // For RAND_bytes(.)

//#define VERBOSE

ZKBPP::ZKBPP() {
  
}

ZKBPP::~ZKBPP() {
  
}

void ZKBPP::init(uint32 party_size, uint32 num_iterations, CircuitContainer* cc, bool print_result) {
  this->party_size_ = party_size;
  this->num_iterations_ = num_iterations;
  this->circuit_ = cc;
  // Get cipher params
  cc->getParams(&(this->circuit_value_size_), &(this->random_tape_size_), &(this->circuit_key_size_), &(this->circuit_gate_size_), &(this->circuit_num_view_gates_));

  // Hardcode this->hash_size_ to SHA256_DIGEST_LENGTH (32 bytes)
  this->hash_size_ = SHA256_DIGEST_LENGTH;
  // Hardcode this->random_tape_size_ to this->circuit_value_size_
  //this->random_tape_size_ = this->circuit_value_size_;
  this->view_size_ = this->circuit_num_view_gates_ * this->circuit_gate_size_;

  this->last_gensign_time_ = 0;
  this->last_sign_time_ = 0;
  this->last_genverify_time_ = 0;
  this->last_verify_time_ = 0;

  // Utils
  this->print_result_ = print_result;

  // Init circuit
  //circuit->init(gate_size);
}

Proof* ZKBPP::sign(uchar* x) {
  auto gensign_start = std::chrono::high_resolution_clock::now();

  // Maybe use alignas(8) or alignas(32) for the arrays below
  // Create buffer for views
  uchar* buffer_all = new uchar[this->num_iterations_ * this->party_size_ * this->view_size_ + // views
    this->num_iterations_ * this->party_size_ * this->random_tape_size_]; // random tapes
  //uchar sign_views_buffer[this->num_iterations_][this->party_size_][this->view_size_];
  uchar* sign_views_buffer = buffer_all;
  memset(sign_views_buffer, 0, this->num_iterations_ * this->party_size_ * this->view_size_);
  
  // Create buffer for random tapes
  //uchar random_tapes_buffer[this->num_iterations_][this->party_size_ * this->circuit_value_size_];
  uchar* random_tapes_buffer = buffer_all + (this->num_iterations_ * this->party_size_ * this->view_size_);
  //RAND_bytes((uchar*)random_tapes_buffer, this->num_iterations_ * this->party_size_ * this->circuit_value_size_);
  BigIntLib::FillRandom(random_tapes_buffer, (this->circuit_)->getKey(), this->num_iterations_ * this->party_size_ * this->random_tape_size_);
  //BigIntLib::Print(random_tapes_buffer, this->num_iterations_ * this->party_size_ * this->circuit_value_size_);

  // Create data
  uint32 k_View_size = this->random_tape_size_ + (this->circuit_num_view_gates_ * this->circuit_gate_size_);
  ContainerSignData* csd = this->createContainerSignData(sign_views_buffer, random_tapes_buffer);
  ContainerCD* ccd = this->createContainerCD(k_View_size, this->party_size_);
  ContainerA* ca = this->createContainerA();
  auto gensign_stop = std::chrono::high_resolution_clock::now();
  this->last_gensign_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(gensign_stop - gensign_start).count();
  auto sign_start = std::chrono::high_resolution_clock::now();
  // For each iteration (iterations are independent from each other):
  //#pragma omp parallel for
  uint32 hash_data_size = this->random_tape_size_ + this->view_size_;
  alignas(32) uchar hash_data[hash_data_size];
  //uchar* hash_data = new uchar[hash_data_size];
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    // Call the circuit, pass SignData and let the circuit fill all Views (and the k's)
    //this->circuit_->evaluateSign(x, x_size, csd->sds_[i]);
    (this->circuit_->runSign)(x, csd->sds_[i]);
    // Add to CD and A
    this->fillCDSign(ccd, i, csd->sds_[i], hash_data);
    this->fillASign(ca, i, csd->sds_[i], ccd);
  }
  //delete[] hash_data;
  //auto sign_stop = std::chrono::high_resolution_clock::now();

  // Commitment
  this->commitment_ = ccd;
  // Allocate memory for the proof p, containing one b_i and one z_i for each iteration
  Proof* p = this->createProof();
  // Create the callenge E, being the hash value of the concatenation of all a_i's
  this->buildChallengeHash(p->e_, ca);
  // For each iteration, the hash value E should create an e_i in {0, 1, 2}
  // For each iteration (iterations are independent from each other):
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    // 1. Create b_i = (y_e+2, C_e+2), where y and C are this iteration's values and e is in {0, 1, 2} according to the challenge (hash value) E
    // 2. Create z_i for this iteration and use this iteration's View, k values and (if e = 2) x_3
    // Remarks: SignData* from this iteration still contains everything needed for z_i (all views, all random tapes and x_3) so according to e_i, move correct pointers into new z_i and delete others from SignData's memory (e.g. two unneeded views for each iteration)
    this->fillProof(p, i, csd->sds_[i], ccd);
  }

  auto sign_stop = std::chrono::high_resolution_clock::now();
  this->last_sign_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(sign_stop - sign_start).count();

  // Delete unneeded resources
  this->destroyContainerSignData(csd);
  this->destroyContainerA(ca);

  delete[] buffer_all;
  return p;
}

bool ZKBPP::verify(Proof* p, uchar* x, uchar* y) {

  //this->printProof(p, false);
  auto genverify_start = std::chrono::high_resolution_clock::now();

  // Create buffer for views
  uchar verify_views_buffer[this->num_iterations_][this->view_size_];
  memset(verify_views_buffer, 0, this->num_iterations_ * this->view_size_);

  bool ret_val;
  uint32 k_View_size = this->random_tape_size_ + (this->circuit_num_view_gates_ * this->circuit_gate_size_);
  ContainerVerifyData* cvd = this->createContainerVerifyData(verify_views_buffer);
  ContainerCD* ccd_verify = this->createContainerCD(k_View_size, this->party_size_ - 1);
  ContainerA* ca_verify = this->createContainerA();
  auto genverify_stop = std::chrono::high_resolution_clock::now();
  this->last_genverify_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(genverify_stop - genverify_start).count();
  auto verify_start = std::chrono::high_resolution_clock::now();
  // For each iteration (iterations are independent from each other):
  uint32 hash_data_size = this->random_tape_size_ + this->view_size_;
  alignas(32) uchar hash_data[hash_data_size];
  //uchar* hash_data = new uchar[hash_data_size];
  //#pragma omp parallel for
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    // Call the circuit, pass VerifyData and let the circuit fill the chosen View
    //this->circuit_->evaluateVerify(p, y, cvd->vds_[i], i); // Maybe better to pass only relevant z_i and not the entire p
    (this->circuit_->runVerify)(p, x, y, cvd->vds_[i], i);
    // Add to CD and A
    this->fillCDVerify(ccd_verify, i, p, cvd->vds_[i], hash_data);
    this->fillAVerify(ca_verify, i, p, cvd->vds_[i], ccd_verify);
  }
  //delete[] hash_data;

  // Create the callenge E, compare with challenge contained in proof
  uchar* challenge_prime = new uchar[this->hash_size_];
  this->buildChallengeHash(challenge_prime, ca_verify);
  if(this->print_result_ == true) {
    std::cout << "[ZKBPP] challenge: " << std::endl;
    this->printDataAsHex(p->e_, this->hash_size_, true);
    std::cout << "[ZKBPP] challenge': " << std::endl;
    this->printDataAsHex(challenge_prime, this->hash_size_, true);
  }
  uint32 result = memcmp(p->e_, challenge_prime, this->hash_size_);
  if(result == 0)
    ret_val = true;
  else
    ret_val = false;

  auto verify_stop = std::chrono::high_resolution_clock::now();
  this->last_verify_time_ = std::chrono::duration_cast<std::chrono::nanoseconds>(verify_stop - verify_start).count();

  // Clean up
  this->destroyContainerVerifyData(cvd);
  this->destroyContainerCD(ccd_verify, this->party_size_ - 1);
  this->destroyContainerA(ca_verify);
  delete[] challenge_prime;
  this->destroyContainerCD(this->commitment_, this->party_size_);
  this->destroyProof(p);

  return ret_val;
}

uint64 ZKBPP::getLastGenSignNS() {
  return this->last_gensign_time_;
}

uint32 ZKBPP::getLastSignNS() {
  return this->last_sign_time_;
}

uint64 ZKBPP::getLastGenVerifyNS() {
  return this->last_genverify_time_;
}

uint32 ZKBPP::getLastVerifyNS() {
  return this->last_verify_time_;
}

uint32 ZKBPP::getExpectedSignatureSizeBits() {
  uint32 hash_bits = this->hash_size_ * 8;
  uint32 value_bits = this->circuit_->getCipherNumBranches() * BigIntLib::field_size_bits_;
  return (this->num_iterations_ * (hash_bits + 2 * value_bits + log2(3) + (2.0/3.0) * value_bits + this->circuit_->getViewSizeBits()));
}

View* ZKBPP::createView() {
  View* v = new View;
  v->values_ = new uchar[this->view_size_];
  memset(v->values_, 0, this->view_size_);
  v->num_gates_ = this->circuit_num_view_gates_;
  v->gate_size_ = this->circuit_gate_size_;
  return v;
}

void ZKBPP::destroyView(View* v) {
  /*
  uint32 num_gates = v->num_gates_;
  for(uint32 i = 0; i < num_gates; i++) {
    delete[] v->values_[i];
  }
  */
  delete[] v->values_;
  delete v;
}

SignData* ZKBPP::createSignData(void* sign_views_buffer, void* random_tapes_buffer, uint32 iteration) {
  uchar (*pointer_1)[this->party_size_][this->view_size_] = (uchar (*)[this->party_size_][this->view_size_]) sign_views_buffer;
  uchar (*pointer_2)[this->party_size_ * this->random_tape_size_] = (uchar (*)[this->party_size_ * this->random_tape_size_]) random_tapes_buffer;
  SignData* sign_data = new SignData;
  sign_data->x_3_ = new uchar[this->circuit_value_size_];
  sign_data->y_shares_ = new uchar[this->party_size_ * this->circuit_value_size_];
  //sign_data->random_tapes_ = new uchar*[this->party_size_];
  sign_data->random_tapes_ = pointer_2[iteration];
  sign_data->random_tapes_hashs_ = new uchar[(this->party_size_ - 1) * this->hash_size_];
  SHA256Dash(sign_data->random_tapes_hashs_, sign_data->random_tapes_, this->random_tape_size_);
  SHA256Dash(sign_data->random_tapes_hashs_ + this->hash_size_, sign_data->random_tapes_ + this->random_tape_size_, this->random_tape_size_);
  
  sign_data->views_ = new uchar*[this->party_size_];
  sign_data->views_[0] = pointer_1[iteration][0];
  sign_data->views_[1] = pointer_1[iteration][1];
  sign_data->views_[2] = pointer_1[iteration][2];
  
  sign_data->y_ = new uchar[this->circuit_value_size_];
  return sign_data;
}

void ZKBPP::destroySignData(SignData* sign_data) {
  delete[] sign_data->x_3_;
  delete[] sign_data->y_shares_;
  //delete[] sign_data->random_tapes_;
  delete[] sign_data->random_tapes_hashs_;
  delete[] sign_data->views_;
  delete[] sign_data->y_;
  delete sign_data;
}

VerifyData* ZKBPP::createVerifyData(void* verify_views_buffer, uint32 iteration) {
  uchar (*pointer_1)[this->view_size_] = (uchar (*)[this->view_size_]) verify_views_buffer;
  VerifyData* verify_data = new VerifyData;
  //verify_data->view_ = new uchar[this->view_size_];
  //memset(verify_data->view_, 0, this->view_size_);
  verify_data->view_ = pointer_1[iteration];
  verify_data->y_share_ = new uchar[this->circuit_value_size_];
  verify_data->y_e2_ = new uchar[this->circuit_value_size_];
  return verify_data;
}

void ZKBPP::destroyVerifyData(VerifyData* verify_data) {
  //delete[] verify_data->view_;
  delete[] verify_data->y_share_;
  delete[] verify_data->y_e2_;
  delete verify_data;
}

ContainerSignData* ZKBPP::createContainerSignData(void* sign_views_buffer, void* random_tapes_buffer) {
  ContainerSignData* csd = new ContainerSignData;
  csd->sds_ = new SignData*[this->num_iterations_];
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    csd->sds_[i] = this->createSignData(sign_views_buffer, random_tapes_buffer, i);
  }
  return csd;
}

void ZKBPP::destroyContainerSignData(ContainerSignData* csd) {
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    destroySignData(csd->sds_[i]);
  }
  delete[] csd->sds_;
  delete csd;
}

ContainerVerifyData* ZKBPP::createContainerVerifyData(void* verify_views_buffer) {
  ContainerVerifyData* cvd = new ContainerVerifyData;
  cvd->vds_ = new VerifyData*[this->num_iterations_];
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    cvd->vds_[i] = this->createVerifyData(verify_views_buffer, i);
  }
  return cvd;
}

void ZKBPP::destroyContainerVerifyData(ContainerVerifyData* cvd) {
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    destroyVerifyData(cvd->vds_[i]);
  }
  delete[] cvd->vds_;
  delete cvd;
}

C* ZKBPP::createC() {
  C* c = new C;
  c->H_k_View_ = new uchar[this->hash_size_];
  return c;
}

void ZKBPP::destroyC(C* c) {
  delete[] c->H_k_View_;
  delete c;
}

D* ZKBPP::createD(uint32 k_View_size) {
  D* d = new D;
  d->k_View_ = new uchar[k_View_size];
  return d;
}

void ZKBPP::destroyD(D* d) {
  delete[] d->k_View_;
  delete d;
}

ContainerCD* ZKBPP::createContainerCD(uint32 k_View_size, uint32 party_size) {
  ContainerCD* ccd = new ContainerCD;
  ccd->Cs_ = new C**[party_size];
  ccd->Ds_ = new D**[party_size];
  // (party_size * iterations) entries
  for(uint32 i = 0; i < party_size; i++) {
    ccd->Cs_[i] = new C*[this->num_iterations_];
    ccd->Ds_[i] = new D*[this->num_iterations_];
    for(uint32 j = 0; j < this->num_iterations_; j++) {
      ccd->Cs_[i][j] = this->createC();
      ccd->Ds_[i][j] = this->createD(k_View_size);
    }
  }
  return ccd;
}

void ZKBPP::destroyContainerCD(ContainerCD* ccd, uint32 party_size) {
  // (party_size * iterations) entries
  for(uint32 i = 0; i < party_size; i++) {
    for(uint32 j = 0; j < this->num_iterations_; j++) {
      this->destroyC(ccd->Cs_[i][j]);
      this->destroyD(ccd->Ds_[i][j]);
    }
    delete[] ccd->Cs_[i];
    delete[] ccd->Ds_[i];
  }
  delete[] ccd->Cs_;
  delete[] ccd->Ds_;
  delete ccd;
}

A* ZKBPP::createA() {
  A* a = new A;
  a->ys_C_hashs_ = new uchar[(this->party_size_ * this->circuit_value_size_) + (this->party_size_ * this->hash_size_)];
  return a;
}

void ZKBPP::destroyA(A* a) {
  delete[] a->ys_C_hashs_;
  delete a;
}

ContainerA* ZKBPP::createContainerA() {
  ContainerA* ca = new ContainerA;
  ca->as_ = new A*[this->num_iterations_];
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    ca->as_[i] = this->createA();
  }
  return ca;
}

void ZKBPP::destroyContainerA(ContainerA* ca) {
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    this->destroyA(ca->as_[i]);
  }
  delete[] ca->as_;
  delete ca;
}

B* ZKBPP::createB() {
  B* b = new B;
  b->y_e2_ = new uchar[this->circuit_value_size_];
  b->H_k_View_ = new uchar[this->hash_size_];
  return b;
}

void ZKBPP::destroyB(B* b) {
  delete[] b->y_e2_;
  delete[] b->H_k_View_;
  delete b;
}

Z* ZKBPP::createZ(bool create_x_3) {
  Z* z = new Z;
  z->view_ = new uchar[this->view_size_];
  z->k_1_ = new uchar[this->random_tape_size_];
  z->k_2_ = new uchar[this->random_tape_size_];
  z->k_1_hash_ = NULL;
  z->k_2_hash_ = NULL;
  if(create_x_3 == true)
    z->x_3_ = new uchar[this->circuit_value_size_];
  else
    z->x_3_ = NULL;
  z->key_shares_ = new uchar*[this->party_size_ - 1];
  for(uint32 i = 0; i < this->party_size_ - 1; i++) {
    z->key_shares_[i] = new uchar[this->circuit_key_size_];
  }
  z->y_share_ = new uchar[this->circuit_value_size_];
  return z;
}

void ZKBPP::destroyZ(Z* z) {
  delete[] z->view_;
  delete[] z->k_1_;
  delete[] z->k_2_;
  if(z->k_1_hash_ != NULL) delete[] z->k_1_hash_;
  if(z->k_2_hash_ != NULL) delete[] z->k_2_hash_;
  if(z->x_3_ != NULL)
    delete[] z->x_3_;
  for(uint32 i = 0; i < this->party_size_ - 1; i++) {
    delete[] z->key_shares_[i];
  }
  delete[] z->key_shares_;
  delete[] z->y_share_;
  delete z;
}

Proof* ZKBPP::createProof() {
  Proof* p = new Proof;
  p->num_iterations_ = this->num_iterations_;
  p->e_ = new uchar[this->hash_size_];
  //p->bs_ = new B*[this->num_iterations_];
  p->y_e2_ = new uchar*[this->num_iterations_];
  p->H_k_View_ = new uchar*[this->num_iterations_];
  p->zs_ = new Z*[this->num_iterations_];
  p->ys_ = new uchar*[this->num_iterations_];
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    //p->bs_[i] = this->createB();
    p->y_e2_[i] = new uchar[this->circuit_value_size_];
    p->H_k_View_[i] = new uchar[this->hash_size_];
    p->zs_[i] = NULL; // Create dynamically later, because x_3 is not always needed!
    p->ys_[i] = new uchar[this->circuit_value_size_];
  }
  return p;
}

void ZKBPP::destroyProof(Proof* p) {
  delete[] p->e_;
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    //this->destroyB(p->bs_[i]);
    delete[] p->y_e2_[i];
    delete[] p->H_k_View_[i];
    this->destroyZ(p->zs_[i]);
    delete[] p->ys_[i];
  }
  delete[] p->y_e2_;
  delete[] p->H_k_View_;
  delete[] p->zs_;
  delete[] p->ys_;
  delete p;
}

void ZKBPP::fillCDSign(ContainerCD* ccd, uint32 iteration, SignData* sign_data, uchar* hash_data) {
  // All <party_size> needed views and tapes are in sign_data
  uint32 hash_data_size = this->random_tape_size_ + this->view_size_;
  //uchar* hash_data = new uchar[hash_data_size];
  //uchar hash_data[hash_data_size] __attribute__ ((aligned (32)));
  for(uint32 i = 0; i < this->party_size_; i++) {
    //uint32 offset = 0;
    memcpy(hash_data, sign_data->random_tapes_ + (i * this->random_tape_size_), this->random_tape_size_);
    memcpy(hash_data + this->random_tape_size_, sign_data->views_[i], this->view_size_);
    //memcpy(hash_data + this->random_tape_size_, sign_data->views_ + (i * this->view_size_), this->view_size_);
    // Build hash and store in iteration's C
    this->SHA256Prime((ccd->Cs_[i][iteration])->H_k_View_, hash_data, hash_data_size);
    // Store concatenation in iteration's D
    memcpy((ccd->Ds_[i][iteration])->k_View_, hash_data, hash_data_size);
  }

  //delete[] hash_data;
}

void ZKBPP::fillCDVerify(ContainerCD* ccd, uint32 iteration, Proof* p, VerifyData* verify_data, uchar* hash_data) {
  // All <party_size - 1> needed views and tapes are in p and verify_data
  uchar* views_temp[this->party_size_ - 1];
  uchar* r_tapes_temp[this->party_size_ - 1];
  views_temp[0] = verify_data->view_; // Calculated by Circuit->evaluateVerify(.)
  views_temp[1] = (p->zs_[iteration])->view_;
  r_tapes_temp[0] = (p->zs_[iteration])->k_1_;
  r_tapes_temp[1] = (p->zs_[iteration])->k_2_;
  uint32 hash_data_size = this->random_tape_size_ + this->view_size_;
  //uchar* hash_data = new uchar[hash_data_size];
  //uchar hash_data[hash_data_size];
  for(uint32 i = 0; i < (this->party_size_ - 1); i++) {
    //uint32 offset = 0;
    memcpy(hash_data, r_tapes_temp[i], this->random_tape_size_);
    memcpy(hash_data + this->random_tape_size_, views_temp[i], this->view_size_);
    // Build hash and store in iteration's C
    this->SHA256Prime((ccd->Cs_[i][iteration])->H_k_View_, hash_data, hash_data_size);
    // Store concatenation in iteration's D
    memcpy((ccd->Ds_[i][iteration])->k_View_, hash_data, hash_data_size);
  }

  //delete[] hash_data;
}

// fills a_i = [y_1, C_1, y_2, C_2, y_3, C_3]_i
void ZKBPP::fillASign(ContainerA* ca, uint32 iteration, SignData* sign_data, ContainerCD* ccd) {
  uint32 offset = 0;
  for(uint32 i = 0; i < this->party_size_; i++) {
    memcpy(ca->as_[iteration]->ys_C_hashs_ + offset, sign_data->y_shares_ + (i * this->circuit_value_size_), this->circuit_value_size_);
    offset += this->circuit_value_size_;
    memcpy(ca->as_[iteration]->ys_C_hashs_ + offset, (ccd->Cs_[i][iteration])->H_k_View_, this->hash_size_); // C value
    offset += this->hash_size_;
  }
}

// fills a_i = [y_1, C_1, y_2, C_2, y_3, C_3]_i
void ZKBPP::fillAVerify(ContainerA* ca, uint32 iteration, Proof* p, VerifyData* verify_data, ContainerCD* ccd) {
  uint32 e = p->e_[iteration % this->hash_size_] % this->party_size_; // TODO
  
  uchar* ys_temp[this->party_size_];
  uchar* C_hashs_temp[this->party_size_];
  ys_temp[e] = verify_data->y_share_; // Last output from calculated View
  ys_temp[(e + 1) % this->party_size_] = (p->zs_[iteration])->y_share_; // Last output from View contained in Proof
  ys_temp[(e + 2) % this->party_size_] = verify_data->y_e2_; // Part of VerifyData (TODO or use y_e+2 from b_i in proof?)
  C_hashs_temp[e] = ccd->Cs_[0][iteration]->H_k_View_; // Calculated C hash value from calculated View
  C_hashs_temp[(e + 1) % this->party_size_] = ccd->Cs_[1][iteration]->H_k_View_; // Calculated C hash value from View contained in Proof
  C_hashs_temp[(e + 2) % this->party_size_] = p->H_k_View_[iteration];
  uint32 offset = 0;
  for(uint32 i = 0; i < this->party_size_; i++) {
    memcpy(ca->as_[iteration]->ys_C_hashs_ + offset, ys_temp[i], this->circuit_value_size_); // y value
    offset += this->circuit_value_size_;
    memcpy(ca->as_[iteration]->ys_C_hashs_ + offset, C_hashs_temp[i], this->hash_size_); // C value
    offset += this->hash_size_;
  }
}

void ZKBPP::fillProof(Proof* p, uint32 iteration, SignData* sign_data, ContainerCD* ccd) {
  uint32 e = p->e_[iteration % this->hash_size_] % this->party_size_; // TODO

  // Add y_i (should all be the same)
  memcpy(p->ys_[iteration], sign_data->y_, this->circuit_value_size_);

  // Add b_i
  memcpy(p->y_e2_[iteration], sign_data->y_shares_ + (((e + 2) % this->party_size_) * this->circuit_value_size_), this->circuit_value_size_); // y_e+2
  memcpy(p->H_k_View_[iteration], ccd->Cs_[(e + 2) % this->party_size_][iteration]->H_k_View_, this->hash_size_); // C_e+2

  // Create and add z_i
  Z* z;
  if(e == 0) {
    p->zs_[iteration] = this->createZ(false);
  }
  else {
    p->zs_[iteration] = this->createZ(true);
    memcpy((p->zs_[iteration])->x_3_, sign_data->x_3_, this->circuit_value_size_); // x_3
  }
  
  // REMARK: If e in {1, 2}, only one hash value is actually needed! This can be further optimized.
  memcpy((p->zs_[iteration])->k_1_, sign_data->random_tapes_ + (e * this->random_tape_size_), this->random_tape_size_); // k_e
  memcpy((p->zs_[iteration])->k_2_, sign_data->random_tapes_ + (((e + 1) % this->party_size_) * this->random_tape_size_), this->random_tape_size_); // k_e+1

  // Needed hash values
  if(e == 0) {
    (p->zs_[iteration])->k_1_hash_ = new uchar[this->hash_size_];
    (p->zs_[iteration])->k_2_hash_ = new uchar[this->hash_size_];
    memcpy((p->zs_[iteration])->k_1_hash_, sign_data->random_tapes_hashs_, this->hash_size_);
    memcpy((p->zs_[iteration])->k_2_hash_, sign_data->random_tapes_hashs_ + this->hash_size_, this->hash_size_);
  }
  else if(e == 1) {
    (p->zs_[iteration])->k_2_hash_ = new uchar[this->hash_size_];
    memcpy((p->zs_[iteration])->k_2_hash_, sign_data->random_tapes_hashs_ + this->hash_size_, this->hash_size_);
  }
  else { // e == 2
    (p->zs_[iteration])->k_1_hash_ = new uchar[this->hash_size_];
    memcpy((p->zs_[iteration])->k_1_hash_, sign_data->random_tapes_hashs_, this->hash_size_);
  }

  // y share
  memcpy((p->zs_[iteration])->y_share_, sign_data->y_shares_ + (((e + 1) % this->party_size_) * this->circuit_value_size_), this->circuit_value_size_);

  memcpy((p->zs_[iteration])->view_, sign_data->views_[(e + 1) % this->party_size_], this->view_size_);
  //memcpy((p->zs_[iteration])->view_, sign_data->views_ + (((e + 1) % this->party_size_) * this->view_size_), this->view_size_);
  //sign_data->views_[(e + 1) % this->party_size_] = NULL; // Set to null, so it's not destroyed by destroySignData, but later with destroyProof
}

void ZKBPP::buildChallengeHash(uchar* destination, ContainerA* ca) {
  uint32 a_size = this->party_size_ * (this->circuit_value_size_ + this->hash_size_);
  uint32 hash_data_size = this->num_iterations_ * a_size; // This is a lot...
  uchar* hash_data = new uchar[hash_data_size];
  uint32 offset = 0;
  for(uint32 i = 0; i < this->num_iterations_; i++) {
    memcpy(hash_data + offset, ca->as_[i]->ys_C_hashs_, a_size);
    offset += a_size;
  }
  // Build challenge hash and add to proof
  SHA256(destination, hash_data, hash_data_size);
  delete[] hash_data;
}

void ZKBPP::SHA256(uchar* destination, uchar* data, uint32 data_size) {
  SHA256_CTX sha256;
  SHA256_Init(&sha256);
  SHA256_Update(&sha256, data, data_size);
  SHA256_Final(destination, &sha256);
}

void ZKBPP::SHA256Prime(uchar* destination, uchar* data, uint32 data_size) { // TODO
  data[data_size - 1] += 1;
  this->SHA256(destination, data, data_size);
  data[data_size - 1] -= 1;
}

void ZKBPP::SHA256Dash(uchar* destination, uchar* data, uint32 data_size) { // TODO
  data[data_size - 1] += 2;
  this->SHA256(destination, data, data_size);
  data[data_size - 1] -= 2;
}

void ZKBPP::printView(uchar* view) {
  for(uint32 j = 0; j < this->circuit_num_view_gates_; j++) {
    std::cout << "Gate " << j << ": ";
    this->printDataAsHex(view + (j * this->circuit_gate_size_), this->circuit_gate_size_, true);
  }
}

void ZKBPP::printProof(Proof* p, bool print_view) {
  bool format = true;
  std::cout << "--- Proof ---" << std::endl;
  std::cout << "Challenge:" << std::endl;
  this->printDataAsHex(p->e_, this->hash_size_, format);
  for(uint32 i = 0; i < p->num_iterations_; i++) {
    std::cout << "--- Iteration " << i << " ---" << std::endl;
    std::cout << "y_i (this is only for testing purposes, should be the same in each iteration):" << std::endl;
    this->printDataAsHex(p->ys_[i], this->circuit_value_size_, format);
    std::cout << "(b) y_e+2:" << std::endl;
    this->printDataAsHex(p->y_e2_[i], this->circuit_value_size_, format);
    std::cout << "(b) H'(k, View):" << std::endl;
    this->printDataAsHex(p->H_k_View_[i], this->hash_size_, format);
    std::cout << "(z) k_e:" << std::endl;
    this->printDataAsHex(p->zs_[i]->k_1_, this->random_tape_size_, format);
    std::cout << "(z) k_e+1:" << std::endl;
    this->printDataAsHex(p->zs_[i]->k_2_, this->random_tape_size_, format);
    std::cout << "(z) x_3:" << std::endl;
    this->printDataAsHex(p->zs_[i]->x_3_, this->circuit_value_size_, format);
    if(print_view == true) {
      std::cout << "(z) View_e+1" << std::endl;
      this->printView(p->zs_[i]->view_);
    }
  }
  std::cout << "-------------" << std::endl;
}

void ZKBPP::printDataAsHex(uchar* data, uint32 data_size, bool format) {
  for(uint32 i = 0; i < data_size; i++) {
    std::cout << std::setfill('0') << std::setw(2) << std::hex << (uint32)data[i];
    if(format == true)
      std::cout << (((i + 1) % 32 == 0) ? "\n" : " ");
  }
  std::cout << std::dec << std::endl;
}
