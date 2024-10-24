#ifndef __FITNESS_LANDSCAPES
#define __FITNESS_LANDSCAPES
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

// there's a way to do this with #define _USE_MATH_DEFINES but it isn't worth debugging it
#define M_PI 3.14159265358979323846

float boxmuller() {
    float u1 = (rand() + 0.0)/RAND_MAX;
    float u2 = (rand() + 0.0)/RAND_MAX;
    return sqrt(-2*log(u2))*cos(2*M_PI*u1);
}

typedef struct node {
    uint32_t genotype; // bitstring
    uint32_t rank;     // integer, least fit is rank 1
    float fitness;     // real number, higher is more fit
} node;

typedef struct K_neighborhood {
    uint32_t* B; // K-length list of the loci in the order they are used to generate the sub-genotype
    float* fitnesses; // 2^K-length list of the fitness values assigned to each sub-genotype
} K_neighborhood;

typedef enum B_styles {BLOCKED, ADJACENT, RANDOM} B_styles;
typedef enum K_styles {CLASSIC, RMF_1, RMF_2, RMF_3} K_styles;

typedef struct NK_landscape {
    uint32_t N; // number of genes
    uint32_t K; // interaction number
    B_styles bst; // method of neighborhood forming
    K_styles kst; // fitness function shape
    node* nodes; // 2^N-length list of nodes
    K_neighborhood* nbs; // N-length list of neighborhoods
} NK_landscape;

typedef struct HNK_landscape {
    uint32_t H; // number of control genes
    uint32_t N; // number of normal genes
    uint32_t K; // interaction number
    B_styles bst;
    node* nodes; // H+N-length bitstring nodes
    uint32_t* G; // control interaction sets (expressed as N length bit-strings)
    NK_landscape gamma_landscape; // "masked" genotype landscape
} HNK_landscape;

uint32_t hammd(uint32_t s1, uint32_t s2) {
    uint32_t diff = s1 ^ s2;
    uint32_t c = 0;
    while (diff > 0) {
        c += (diff & 1);
        diff >>= 1;
    }
    return c;
}
#endif