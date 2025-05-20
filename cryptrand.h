#ifndef __CRYPTRAND__H
#define __CRYPTRAND__H

#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <stdint.h>
#include <functional>
#include <string>

// Include OpenSSL for SHA256 (you'll need to install it: e.g., libssl-dev on Debian/Ubuntu)
// #include <openssl/sha.h>

// # gcc -O3 gen_random.c -o gen_random -lcrypto

uint64_t generate_random_64bit_from_timestamp();

#endif