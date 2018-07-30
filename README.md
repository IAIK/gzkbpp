# Generalized ZKB++
This is an implementation of the [ZKB++ proof system](https://eprint.iacr.org/2017/279.pdf), where the inner circuit can easily be exchanged.

## Dependencies
The OpenSSL library is required.

## Compilation
Use the provided `_compile_zkbpp.sh` file to compile the program:
1. `chmod +x _compile_zkbpp.sh`
2. `./_compile_zkbpp.sh`

The implementation was tested on CentOS 7.5 and by using Valgrind 3.13.0.

## Usage
After compilation, the program can be used with the following parameters:

`./zkboo_test <field_size> <num_branches> <field_type> <cipher_type> <print_result>`

where

* `field_size` specifies the field size in bits,
* `num_branches` specifies the number of branches,
* `field_type` specifies the type of the field (0 for prime field, 1 for binary field),
* `cipher_type` specifies the cipher being used, and
* `print_result` specifies whether results should be printed to the console.


## Circuits
### Adding Circuits
New circuits can be implemented by adding three specific methods to the `CircuitContainer` class (exchange "Instance" with the name of the circuit):
1. `initInstance(.)`

    This method initializes a new circuit. Typical values for block ciphers include the number of rounds or the number of branches for Feistel-based constructions. Any necessary precompution (e.g. round keys, round constants) should also be done in this method. Note that in ZKB++, each multiplication gate requires a set of random values, and therefore the number of multiplication gates must also be set here. Moreover, the necessary function pointers are also set in this method.
2. `directInstance(.)`

    This method is a direct evaluation of the specified circuit, without any multi-party computations, shares or random numbers for the multiplications.
3. `circuitInstance(.)`

    The MPC version of the circuit is implemented in this method. The function pointers previously defined in `beforeSign(.)` and `beforeVerify(.)` have to be used here, because this method is called both during signature generation and signature verification.

### Circuit Gates
A few circuit gates are already implemented. Common circuit gates include addition and multiplication, both with shared or constant values. Note that some circuit gates are implemented differently for signature generation and signature verification. If additional gates are needed, they can be added easily and by using the same techniques found in e.g. `addSharedSign(.)` or `addSharedVerify(.)`.

### Field Arithmetic
The `BigIntLib` class includes methods for computations in a set of predefined finite fields, in particular prime fields and binary fields. New methods can easily be added for finite fields of different sizes. Currently, the Solinas reduction [1] is used in prime fields, and a fast word-wise reduction [2] is used in binary fields.

As an example circuit, [MiMC](https://eprint.iacr.org/2016/492.pdf) is included.

[1] Jerome A. Solinas: Generalized Mersenne Numbers  
[2] Darrel Hankerson, Alfred Menezes and Scott Vanstone: Guide to Elliptic Curve Cryptography
