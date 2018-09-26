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

#include <iostream>
#include "ZKBPP.h"
#include "CircuitContainer.h"
#include <cstring>
#include <cmath> // ceil(.)
#include <openssl/rand.h> // For RAND_bytes(.)
#include <stdio.h>

// Creates plaintext value
void write_plaintext(uchar* plaintext, uint32 plaintext_size) {
  memset(plaintext, 0, plaintext_size);
  for(uint32 i = 0; i < plaintext_size; i++) {
    plaintext[i] = (uchar)(i); // Write some values
  }
}

int main(int argc, char** argv)
{
  std::cout << "Starting" << std::endl;

  // No error handling here...
  if(argc != 2 && argc != 6) {
    std::cout << "Usage: ./program <num_field_bits> <num_branches> <field_type> <cipher_type> <zkbpp_print_result>" << std::endl;
    std::cout << "Enter ./program help for possible choices." << std::endl;
    return 1;
  }

  std::string param(argv[1]);
  if(param == "help") {
    std::cout << "Example configurations: " << std::endl;
    std::cout << "./program 256 1 0 1 0: Uses MiMC-256 (prime field, print_result = false)" << std::endl;
    std::cout << "./program 272 1 0 1 1: Uses MiMC-272 (prime field, print_result = true)" << std::endl;
    return 0;
  }

  cpu_set_t set;
  CPU_ZERO(&set);        // clear cpu mask
  CPU_SET(0, &set);      // set cpu 0
  sched_setaffinity(0, sizeof(cpu_set_t), &set);  // 0 is the calling process

  uint32 field_bits = atoi(argv[1]);
  uint32 num_branches = atoi(argv[2]);
  uint32 value_size = ceil((field_bits * num_branches) / 8.0); // For value sizes mod 8 != 0: use uint32 value_size = ceil((branch_bits * num_branches) / 8.0);
  uint32 random_tape_size = 16;
  uint32 key_size = value_size;
  uint32 field_type = (bool)(atoi(argv[3]));
  uint32 cipher_type = atoi(argv[4]);
  bool print_result = (bool)(atoi(argv[5]));

  uchar plaintext[value_size];
  uchar y[value_size];
  RAND_bytes(plaintext, value_size);
  //write_plaintext(plaintext, value_size); // DETERMINISTIC FOR TESTING!
  uint32 party_size = 3;
  uint32 zkbpp_iterations = 438; // This should be 438 later
  
  
  //MiMCp* c = new MiMCp(party_size);
  CircuitContainer* c = new CircuitContainer();
  c->init(value_size, random_tape_size, key_size, field_bits, num_branches, field_type, party_size);
  c->initCipher(cipher_type);
  std::string instance_name = "(" + std::to_string(field_bits) + ", " + std::to_string(num_branches) + ", " + std::to_string(c->getCipherNumRounds()) + ")";
  std::string instance_name_sage = c->getCipherName() + "-" + instance_name;

  c->directEncryption(plaintext, y);

  ZKBPP* zkbpp = new ZKBPP();
  zkbpp->init(party_size, zkbpp_iterations, c, print_result);

  uint32 num_runs = 50;
  uint64 total_gensign = 0;
  uint64 total_sign = 0;
  uint64 total_genverify = 0;
  uint64 total_verify = 0;
  uint64 total_circuit_sign = 0;
  uint64 total_circuit_verify = 0;
  Proof* p;
  bool success;
  for(uint32 i = 0; i < num_runs; i++) {
    //c->randomizeKey();
    //c->directEncryption(plaintext, y);

    p = zkbpp->sign(plaintext);
    success = zkbpp->verify(p, plaintext, y);
    total_gensign += zkbpp->getLastGenSignNS();
    total_sign += zkbpp->getLastSignNS();
    total_genverify += zkbpp->getLastGenVerifyNS();
    total_verify += zkbpp->getLastVerifyNS();
    total_circuit_sign += c->getLastCircuitSignNS(); // Last circuit execution of each run
    total_circuit_verify += c->getLastCircuitVerifyNS(); // Last circuit execution of each run
  }

  // Testing direct function
  uint64 min_cycles = 0xFFFFFFFF;
  uint64 total_direct_call = 0;
  for(uint32 i = 0; i < num_runs; i++) {
    c->randomizeKey();
    c->directEncryption(plaintext, y);
    if(c->getLastDirectCallCycles() < min_cycles) min_cycles = c->getLastDirectCallCycles();
    total_direct_call += c->getLastDirectCallNS();
  }

  // Print out information
  float time_gensign = ((float)total_gensign / num_runs / 1000000);
  float time_sign = ((float)total_sign / num_runs / 1000000);
  float time_sign_total =  time_gensign + time_sign;
  float time_genverify = ((float)total_genverify / num_runs / 1000000);
  float time_verify = ((float)total_verify / num_runs / 1000000);
  float time_verify_total = time_genverify + time_verify;
  float time_circuit_sign = ((float)(total_circuit_sign) / num_runs / 1000000);
  float time_circuit_verify = ((float)(total_circuit_verify) / num_runs / 1000000);
  float time_circuits_sign_all = (((float)(total_circuit_sign / num_runs) / 1000000) * zkbpp_iterations);
  float time_circuits_verify_all = (((float)(total_circuit_verify / num_runs) / 1000000) * zkbpp_iterations);
  float time_direct_call = ((float)(total_direct_call) / num_runs / 1000000);
  std::cout << "--- CONFIGURATION ---" << std::endl;
  std::cout << "Number of test runs: " << num_runs << std::endl;
  std::cout << "--- CIRCUIT ---" << std::endl;
  std::cout << "Field type: " << ((field_type == 0) ? "Prime field" : "Binary field") << std::endl;
  std::cout << "Field size (bits): " << field_bits << std::endl;
  std::cout << "Value size (bytes): " << value_size << std::endl;
  std::cout << "Key size (bytes): " << key_size << std::endl;
  std::cout << "Number of branches: " << num_branches << std::endl;
  std::cout << "Number of cipher rounds: " << c->getCipherNumRounds() << std::endl;
  std::cout << "Number of cipher MUL gates: " << c->getCipherNumMulGates() << std::endl;
  std::cout << "[Cipher] Cycles per byte: " << (min_cycles / 32) << std::endl;
  std::cout << "[Cipher] Average time for direct call: " << time_direct_call << " ms" << std::endl;
  std::cout << "--- TIME ---" << std::endl;
  std::cout << "Average time for gensign: " << time_gensign << " ms" << std::endl;
  std::cout << "Average time for sign: " << time_sign << " ms" << std::endl;
  std::cout << "Average time for sign (total): " << time_sign_total << " ms" << std::endl;
  std::cout << "Average time for genverify: " << time_genverify << " ms" << std::endl;
  std::cout << "Average time for verify: " << time_verify << " ms" << std::endl;
  std::cout << "Average time for verify (total): " << time_verify_total << " ms" << std::endl;
  std::cout << "Average time for circuit sign: " << time_circuit_sign << " ms" << std::endl;
  std::cout << "Average time for circuit verify: " << time_circuit_verify << " ms" << std::endl;
  std::cout << "Average time for all circuits sign: " << time_circuits_sign_all << " ms" << std::endl;
  std::cout << "Average time for all circuits verify: " << time_circuits_verify_all << " ms" << std::endl;
  std::cout << "--- ZKB++ ---" << std::endl;
  std::cout << "Number of ZKB++ iterations: " << zkbpp_iterations << std::endl;
  std::cout << "View size: " << (c->getViewSizeBits() / 8) << " bytes (" << c->getViewSizeBits() << " bits)" << std::endl;
  std::cout << "Expected signature size: ~" << (zkbpp->getExpectedSignatureSizeBits() / 8 / 1024) << " KB" << std::endl;
  if(success == true)
    std::cout << "Result: ACCEPT" << std::endl;
  else
    std::cout << "Result: !!!!!!!!!! REJECT !!!!!!!!!!"  << std::endl;

  // Print table for LaTeX (something like: "{{\texttt {Scheme-N-n}}} & {<number> <unit>} & {<number> <unit>} & {<number> <unit>} \\" for each instance (here for this one))
  // | Scheme Definition | GenSign time (ms) | Sign time (ms) | GenVerify Time (ms) | Verify time (ms) | View size (bits) | Signature size (KB) |
  //printf("{{\\texttt {%s}}} & {%.2f ms} & {%.2f ms} & {%d bits} & {\\({\\approx}\\) %d KB} \\\\\n",
  printf("& \\(%s\\) & \\(%.2f\\) ms & \\(%.2f\\) ms & \\(%.2f\\) ms & \\(%.2f\\) ms & \\(%d\\) bits & \\(\\approx %d\\) KB \\\\ \\cline{2-8}\n",
    instance_name.c_str(), time_gensign, time_sign_total, time_genverify, time_verify_total, c->getViewSizeBits(), (zkbpp->getExpectedSignatureSizeBits() / 8 / 1024));

  // Print data for Sage
  printf("[\"%s\", %.2f, %.2f, %d]\n",
    instance_name_sage.c_str(), time_sign_total, time_verify_total, (zkbpp->getExpectedSignatureSizeBits() / 8 / 1024));

  // Clean up
  delete zkbpp;
  delete c;

  return 0;
}
