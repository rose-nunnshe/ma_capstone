#ifndef __MT19937
#define __MT19937

#include <stdio.h>

void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_float32_full(void);
double genrand_float32_notone(void);

#endif