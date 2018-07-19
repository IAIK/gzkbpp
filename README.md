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
After compilation, the program can be called with following parameters:

`./zkboo_test <field_size> <num_branches> <cipher> <field_type> <print_result>`

where

* `field_size` specifies the field size in bits,
* `num_branches` specifies the number of branches,
* `cipher_type` specifies the cipher being used,
* `field_type` specifies the type of the field (0 for prime field, 1 for binary field), and
* `print_result` specifies whether results should be printed to the console.


## Circuits
As an example, [MiMC](https://eprint.iacr.org/2016/492.pdf) is included.
