#include <fitness_landscapes.h>

const char* kst_to_str(K_styles kst) {
    switch (kst) {
        case CLASSIC:
            return "CLASSIC";
        break;
        case RMF_1:
            return "RMF_Theta=1";
        break;
        case RMF_2:
            return "RMF_Theta=2";
        break;
        case RMF_3:
            return "RMF_Theta=3";
        break;
        default:
            return "N/A";
        break;
    }
}

const char* bst_to_str(B_styles bst) {
    switch (bst) {
        case ADJACENT:
            return "ADJACENT";
        break;
        case BLOCKED:
            return "BLOCKED";
        break;
        case RANDOM:
            return "RANDOM";
        break;
        default:
            return "N/A";
        break;
    }
}

void print_genotype_bitstring(uint32_t s, uint32_t N) {
    char ret[200];
    ret[0] = ((s & 1)==0)?'0':'1';
    for (uint32_t i = 1; i < N; i++) {
        ret[i] = ((s & (1<<i))==0)?'0':'1';
    }
    ret[N] = 0;
    printf("%s", ret);
}

// Allocates space for nodes and neighborhoods arrays, assigns N and K
// does NOT build anything (because we don't want to be stuck leaking memory)
void init_NK(NK_landscape* nk, uint32_t N, uint32_t K) {
    // printf("<init_NK>\n");
    // basic assignments
    nk->N = N;
    nk->K = K;
    nk->nodes = (node*)malloc((1 << N)*sizeof(node));
    nk->nbs = (K_neighborhood*)malloc(N*sizeof(K_neighborhood));

    // make neighborhoods
    for (uint32_t n = 0; n < N; n++) {
        K_neighborhood* nb = nk->nbs + n;
        nb->B = (uint32_t*)malloc(K*sizeof(uint32_t));
        nb->fitnesses = (float*)malloc((1<<K)*sizeof(float));
    }
    // printf("</init_NK>\n");
}

// Deallocates an NK landscape
void deinit_NK(NK_landscape* nk) {
    free(nk->nodes);
    for (uint32_t n = 0; n < nk->N; n++) {
        K_neighborhood* nb = nk->nbs + n;
        free(nb->B);
        free(nb->fitnesses);
    }
    free(nk->nbs);
}

// arg1 and arg2 are node pointers
int compare1(const void* arg1, const void* arg2) {
    float diff1 = ((node*)arg1)->fitness - ((node*)arg2)->fitness;
    if (diff1 < 0) return -1;
    if (diff1 > 0) return 1;
    return 0;
}
int compare2(const void* arg1, const void* arg2) {
    return ((node*)arg1)->genotype - ((node*)arg2)->genotype;
}

float NK_calc_fitness(NK_landscape* nk, uint32_t s, uint8_t verbose) {
    uint32_t N = nk->N;
    uint32_t K = nk->K;
    float fit = 0.0;
    if (verbose) { 
        printf("Calculating fitness of ");
        print_genotype_bitstring(s, N);
        printf("\n");
    }
    for (uint32_t ell = 0; ell < N; ell++) {
        uint32_t subgenotype = 0;
        for (uint32_t k = 0; k < K; k++) {
            // Translation of this: the subgenotype for locus ell
            // is assembled from the K indices of s dictated by B_ell
            // in the order they appeared in B
            subgenotype |= ((s >> (nk->nbs)[ell].B[k]) & 1) << k;
        }
        if (verbose) {
            printf("... subgenotype_%d = ", ell);
            print_genotype_bitstring(subgenotype, K);
            printf(" with fitness contribution %f\n", nk->nbs[ell].fitnesses[subgenotype]);
        }

        // now we take the subgenotype and look up the fitness contribution from the table
        // and add that into the overall
        fit += (nk->nbs)[ell].fitnesses[subgenotype];
    }
    
    return fit;
}

// builds neighborhoods, random fitness values, etc.
// DOES NOT allocate anything (N and K already need to be correctly set up)
void reroll_NK(NK_landscape* nk, K_styles kst, B_styles bst, uint8_t do_rank) {
    // printf("<reroll_NK>\n");
    uint32_t K = nk->K;
    uint32_t N = nk->N;
    nk->kst = kst;
    nk->bst = bst;
    
    uint8_t succeeded = 0;
    while (!succeeded) {
        succeeded = 1;
        // build neighborhoods
        for (uint32_t n = 0; n < N; n++) {
            K_neighborhood* nb = nk->nbs + n;

            for (uint32_t i = 0; i < K; i++) {
                switch (bst) {
                    case BLOCKED:
                        (nb->B)[i] = ((n/K)*K+i) % N; // creates "square" blocks of interaction
                    break;
                    case ADJACENT:
                        (nb->B)[i] = (n-K/2+i + N) % N; // K interactions centered on n (with wrapping)
                    break;
                    case RANDOM:
                        if (i == 0) {
                            (nb->B)[i] = n; // must guarantee n interacts with itself
                        }
                        else {
                            // it's possible we will randomly choose an already-included locus, so
                            // we have to reroll the random selection until we get a new one
                            uint32_t candidate = rand() % N;
                            uint32_t j = 0;
                            while (j < i) {
                                if (candidate != (nb->B)[j]) { // no collision at this index
                                    j++;
                                }
                                else {
                                    // collision, have to start over
                                    candidate = rand() % N;
                                    j = 0;
                                }
                            }
                            
                            // finally have a new index to add
                            (nb->B)[i] = candidate;
                        }
                    break;
                    default: // shouldn't happen
                        printf("Yikes! Illegal B style.\n");
                        exit(-1);
                    break;
                }
            }

            if (kst == CLASSIC) {
                for (uint32_t s = 0; s < (1<<K); s++) {
                    (nb->fitnesses)[s] = (rand() + 0.0)/RAND_MAX;
                }
            }
            else {
                // need to choose a peak genotype
                uint32_t delta = rand() % (1<<K);

                float theta = 1.0;
                if (kst == RMF_2) theta = 2.0;
                if (kst == RMF_3) theta = 3.0;

                for (uint32_t s = 0; s < (1 << K); s++) {
                    (nb->fitnesses)[s] = -theta*hammd(s, delta) + boxmuller();
                }
            }
        }

        // iterate through the genotypes building their fitness
        for (uint32_t s = 0; s < (1<<N); s++) {
            node* current = (nk->nodes) + s;
            current->genotype = s;
            current->fitness = NK_calc_fitness(nk, s, 0); // don't want verbose
        }

        // need to take a second pass to build their rank metadata (if wanted)
        if (do_rank) {
            // this is a grisly way to do this, but it is likely to be comparably fast 
            // to simply finding the minimum fitness over and over again (which is O(2^(2N))).
            // this approach is O(2*2^N*log(2^N)) = ~O(N*2^N) as opposed to O(4^N). Then again,
            // for N <= 8, probably none of this matters.

            // first, sort the nodes by fitness so they are in ascending order
            qsort(nk->nodes, (1<<(nk->N)), sizeof(node), compare1);

            // then, assign their ranks based on that order
            for (uint32_t s = 0; s < (1<<N); s++) {
                (nk->nodes)[s].rank = s+1; // least fit should be rank 1, not rank 0
            }
            
            // error check: no one can be tied
            for (uint32_t s = 1; s < (1<<N); s++) {
                if ((nk->nodes)[s].fitness == (nk->nodes)[s-1].fitness) {
                    /*
                    printf("***** Holy God, it's a tie! (%f, genotype ", nk->nodes[s].fitness); 
                    print_genotype_bitstring(nk->nodes[s].genotype, nk->N);
                    printf(" with %f, genotype ", nk->nodes[s-1].fitness); 
                    print_genotype_bitstring(nk->nodes[s-1].genotype, nk->N); 
                    printf(") *****\n"); 
                    
                    // need to interrogate exactly how this happened
                    NK_calc_fitness(nk, nk->nodes[s].genotype, 1);
                    NK_calc_fitness(nk, nk->nodes[s-1].genotype, 1);
                    */
                    // if we didn't get an injective function, we just have to try again.
                    succeeded = 0; 
                    break;
                }
            }

            // then, sort the nodes by genotype so they're back in index order
            qsort(nk->nodes, (1<<(nk->N)), sizeof(node), compare2);
        }
    }
    // so now, all nodes have a correct fitness and correct rank
    // printf("</reroll_NK>\n");
}

// Allocates nodes and space for underlying NK landscape
// Note: intolerant of H=0 -- if you want to simulate that, use the NK model 
//       with K_style = CLASSIC
void init_HNK(HNK_landscape* hnk, uint32_t H, uint32_t N, uint32_t K) {
    // printf("<init_HNK>\n");
    hnk->H = H;
    hnk->N = N;
    hnk->K = K;

    hnk->G = (uint32_t*)malloc(H*sizeof(uint32_t));
    hnk->nodes = (node*)malloc((1<<(H+N))*sizeof(node));
    hnk->gamma_landscape = (NK_landscape*)malloc(sizeof(NK_landscape));
    
    init_NK(hnk->gamma_landscape, N, K);
    // printf("</init_HNK>\n");
}

void deinit_HNK(HNK_landscape* hnk) {
    deinit_NK(hnk->gamma_landscape);
    free(hnk->gamma_landscape);
    free(hnk->nodes);
    free(hnk->G);
}

uint32_t HNK_sigma_to_gamma(uint32_t s, uint32_t H, uint32_t* G) {
  uint32_t gamma = s >> H;
  for (uint32_t h = 0; h < H; h++) {
      uint8_t H_h = (s >> h) & 1; // current control gene value
      uint32_t G_h = G[h];

      // what this does: everything not in G_h is unaffected
      // everything in G_h is unaffected if H_h = 1
      // everything in G_h is forced to 0 if H_h = 0
      gamma &= (~G_h) | (G_h*H_h);
  }
  return gamma;
}

// builds neighborhoods, random fitness values, etc.
// DOES NOT allocate anything (H, N, and K already need to be correctly set up)
void reroll_HNK(HNK_landscape* hnk, B_styles bst, uint8_t do_rank) {
    // printf("<reroll_HNK>\n");
    uint32_t H = hnk->H;
    uint32_t N = hnk->N;
    uint32_t K = hnk->K;
    hnk->bst = bst;

    // create the control relationships
    for (uint32_t h = 0; h < H; h++) {
        uint32_t G_h = 0;
        for (uint32_t i = 0; i < K; i++) {
            uint32_t candidate = 1 << (rand() % N);
            while (G_h & candidate) { // this bit was picked before
                candidate = 1 << (rand() % N);
            }
            G_h |= candidate;
        }
        (hnk->G)[h] = G_h;
    }

    // build the underlying classic NK
    reroll_NK(hnk->gamma_landscape, CLASSIC, bst, do_rank);

    // figure out what correspondence to give based on masking
    for (uint32_t s = 0; s < (1<<(H+N)); s++) {
        (hnk->nodes)[s].genotype = s;

        // construct gamma
        uint32_t gamma = HNK_sigma_to_gamma(s, H, hnk->G);

        // give this genotype the same fitness and rank as its gamma correspondence
        (hnk->nodes)[s].fitness = hnk->gamma_landscape->nodes[gamma].fitness;
        (hnk->nodes)[s].rank = hnk->gamma_landscape->nodes[gamma].rank; // if do_rank = 0 these are just garbage
    }
    // printf("</reroll_HNK>\n");
}

// recurrence helper for HNK_has_access. Only reason these two functions are separate
// is that this one obviously does not clear backtracking
uint8_t HNK_has_access_recur(HNK_landscape* hnk, uint32_t s, uint8_t* backtracking) {
    if (backtracking[s]) return 0;
    backtracking[s] = 1;
    
    if ((hnk->nodes)[s].rank == (1<<(hnk->N))) {
        // highest possible rank, this is the global optimum (or could access it)
        return 1;
    }

    // look through all 1-distance neighbors to s and recurse into those that are NOT backtracking
    uint32_t mask = 1;
    while (mask < (1<<(hnk->H+hnk->N))) {
        uint32_t next = s ^ mask;
        if (((hnk->nodes)[next].rank >= (hnk->nodes)[s].rank)) {
            if (HNK_has_access_recur(hnk, next, backtracking)) return 1;
        }
        mask <<= 1;
    }

    return 0;
}

// determines if a genotype s has access to the global optimum on hnk
// hnk must have rank data for this to work. We don't allocate backtracking
// here because that would excessively call malloc if we wanted to call this
// many times on the same landscape (which we will).
uint8_t HNK_has_access(HNK_landscape* hnk, uint32_t s, uint8_t* backtracking) {
    for (uint32_t i = 0; i < (1<<(hnk->N+hnk->H)); i++) {
        backtracking[i] = 0;
    }

    return HNK_has_access_recur(hnk, s, backtracking);
}

uint8_t NK_has_access_recur(NK_landscape* nk, uint32_t s, uint8_t* backtracking) {
    if ((nk->nodes)[s].rank == (1<<(nk->N))) {
        // highest possible rank, this is the global optimum (or could access it)
        return 1;
    }

    uint8_t ret = 0;

    // look through all 1-distance neighbors to s and recurse into those that are NOT backtracking
    uint32_t mask = 1;
    while (mask < (1<<(nk->N))) {
        uint32_t next = s ^ mask;
        if (!(backtracking[next]) && 
             ((nk->nodes)[next].fitness >= (nk->nodes)[s].fitness)) {
            backtracking[next] = 1;
            ret |= NK_has_access_recur(nk, next, backtracking);
        }
        if (ret) return 1; // access is access, don't waste time checking extra paths/nodes
        mask <<= 1;
    }

    return ret;
}

// see above HNK functions for explanation of what this and NK_has_access_recur are doing
// they are essentially the same, except a few nuances about datatypes
uint8_t NK_has_access(NK_landscape* nk, uint32_t s, uint8_t* backtracking) {
    for (uint32_t i = 0; i < (1<<(nk->N)); i++) {
        backtracking[i] = 0;
    }

    return NK_has_access_recur(nk, s, backtracking);
}

// determines the number of 1-local optima on an input NK landscape
// this may seem slow, but in fact it seems unlikely there is any easily
// comprehensible way of doing this that is any faster.
uint32_t NK_count_local_optima(NK_landscape* nk) {
    uint32_t ret = 0;

    uint32_t num_nodes = (1<<(nk->N));
    for (uint32_t s = 0; s < num_nodes; s++) {
        uint32_t mask = 1;
        uint8_t failed = 0;
        while (mask < num_nodes) {
            if ((nk->nodes)[s].fitness < (nk->nodes)[s ^ mask].fitness) {
                failed = 1;
                break;
            }

            mask <<= 1;
        }
        if (!failed) {
            ret++;
        }
    }

    return ret;
}

// finds global optimum, requires rank data to work
uint32_t NK_find_gopt(NK_landscape* nk) {
    for (uint32_t s = 0; s < (1<<(nk->N)); s++) {
        if (nk->nodes[s].rank == (1<<(nk->N))) return s;
    }
    printf("NK_find_gopt: no global optimum found\n");
    return -1;
}

uint32_t NK_calc_p_1_recur(NK_landscape* nk, uint8_t* backtracking, uint32_t s) {
    if (backtracking[s]) return 0;
    backtracking[s] = 1;

    uint32_t ret = 1; // this one
    uint32_t N = nk->N;

    for (uint32_t m = 1; m < (1<<N); m<<=1) {
        uint32_t next = s ^ m;
        if (nk->nodes[next].rank < nk->nodes[s].rank) { // rank must be decreasing working away from global opt 
            ret += NK_calc_p_1_recur(nk, backtracking, next);
        }
    }

    return ret;
}

// calculates p_1 by working from the global optimum outward
float NK_calc_p_1(NK_landscape* nk, uint8_t* backtracking) {
    for (uint32_t s = 0; s < (1<<(nk->N)); s++) {
        backtracking[s] = 0;
    }

    return (NK_calc_p_1_recur(nk, backtracking, NK_find_gopt(nk))+0.0)/(1<<(nk->N));
}

float NK_alt_p_1(NK_landscape* nk, uint8_t* backtracking) {
    uint32_t count = 0;
    for (uint32_t s = 0; s < (1<<(nk->N)); s++) {
        if (NK_has_access(nk, s, backtracking)) count++;
    }
    
    return (count + 0.0)/(1<<(nk->N));
}

// slower HNK local optimum function that is used to check fast version
uint32_t HNK_slow_count_local_optima(HNK_landscape* hnk) {
    uint32_t ret = 0;

    uint32_t num_nodes = (1<<(hnk->N));
    for (uint32_t s = 0; s < num_nodes; s++) {
        uint32_t expanded_s = (s << (hnk->H)) + (1 << (hnk->H)) - 1; // fills in 1s in all the H-segment loci
        uint32_t mask = 1<<(hnk->H);
        uint8_t failed = 0;
        while (mask < (num_nodes << (hnk->H))) {
            if ((hnk->nodes)[expanded_s].fitness < (hnk->nodes)[expanded_s ^ mask].fitness) {
                failed = 1;
                break;
            }

            mask <<= 1;
        }
        if (!failed) {
            ret++;
        }
    }

    return ret;
}

// fast local optimum counting function that uses the fact the count is the 
// same as in the gamma landscape
uint32_t HNK_count_local_optima(HNK_landscape* hnk, uint8_t* backtracking) {
    return NK_count_local_optima(hnk->gamma_landscape);
}

// checks if the global optimum of the HNK landscape is fully expressed
// (that is, is it only the global optimum if all H-segment loci are 1, or
// are there tied genotypes with at least one H-segment locus being 0)
uint8_t HNK_check_gopt(HNK_landscape* hnk, uint32_t gopt) {
    // requires gopt is fully expressed, so now we knock out each individual H-segment locus
    // to see if any of them form ties. If none do, then the global optimum has to 
    // be fully expressed. Note that we only need to try the H-segment loci one
    // at a time, because if there were a tie with two of them or more disabled,
    // it would have to imply either that it was also a tie with only one disabled,
    // or that the underlying fitness landscape (gamma landscape) is not injective.

    // a way to see this more intuitively is to understand that, since the gamma landscape
    // is injective, the only way the global optimum can tie with a non-fully expressed genotype
    // is for all the loci controlled by the relevant control gene to already be expressing
    // type 0. If that's the case, then having 2 genes turned off has to be a superset of 
    // an equivalent situation where only 1 gene was turned off.
    for (uint32_t h = 0; h < (hnk->H); h++) {
        uint32_t h_mask = 1 << h;
        uint32_t gopt_mod = gopt ^ h_mask;
        if (hnk->nodes[gopt_mod].fitness == hnk->nodes[gopt].fitness) {
            return 0; // found a tie
        }
    }
    return 1;
}

uint32_t HNK_find_fully_expressed_gopt(HNK_landscape* hnk) {
    uint32_t H = hnk->H;
    uint32_t N = hnk->N;
    uint32_t h_mask = (1 << H)-1;
    float max_fit = hnk->nodes[h_mask].fitness;
    uint32_t max_s = hnk->nodes[h_mask].genotype;

    for (uint32_t s = (1<<H) + h_mask; s < (1<<(N+H)); s += (1<<H)) {
        float this_fit = hnk->nodes[s].fitness;
        if (this_fit > max_fit) {
            max_fit = this_fit;
            max_s = s;
        }
    }

    return max_s;
}

uint32_t HNK_calc_p_1_recur(HNK_landscape* hnk, uint8_t* backtracking, uint32_t s) {
    if (backtracking[s]) return 0;
    backtracking[s] = 1;

    uint32_t ret = 0;
    uint32_t H = hnk->H;
    uint32_t N = hnk->N;
    uint32_t h_mask = (1<<H)-1;
    if ((s & h_mask) == h_mask) {
        ret += 1; // this is a fully-expressed genotype
        printf("[REGULAR] genotype "); print_genotype_bitstring(s, hnk->H+hnk->N); printf(" has access\n");
    }

    // look for neighbors
    for (uint32_t m = 1; m < (1<<(H+N)); m<<=1) {
        uint32_t next = s ^ m;
        if (hnk->nodes[next].fitness <= hnk->nodes[s].fitness) { // fitness must be non-increasing when moving away from optimum
            ret += HNK_calc_p_1_recur(hnk, backtracking, next);
        }
    }

    return ret;
}

// calculate p_1 for a given HNK landscape (requires rank data)
float HNK_calc_p_1(HNK_landscape* hnk, uint8_t* backtracking) {
    // set up anti-backtracking helpers
    for (uint32_t i = 0; i < (1<<(hnk->N+hnk->H)); i++) {
        backtracking[i] = 0;
    }

    uint32_t gopt = HNK_find_fully_expressed_gopt(hnk);
    // elaborate from gopt with non-increasing fitness, only counting fully expressed genotypes
    // but being allowed to visit non-fully expressed genotypes during recursion
    return (HNK_calc_p_1_recur(hnk, backtracking, gopt)+0.0)/(1<<(hnk->N));
}

// over nReps HNK landscapes, count how many had a fully expressed global optimum, how many had a 
// non-fully expressed global optimum, and p_1 for each case
void test_HNK(uint32_t H, uint32_t N, uint32_t K, B_styles bst, uint32_t nReps, uint32_t* fexn, uint32_t* nexn, float* p_1F, float* p_1N) {
    HNK_landscape hnk;
    init_HNK(&hnk, H, N, K);
    uint8_t* backtracking = (uint8_t*)malloc((1<<(H+N))*sizeof(uint8_t));
    (*fexn) = 0;
    (*nexn) = 0;
    (*p_1F) = 0.0;
    (*p_1N) = 0.0;
    for (uint32_t i = 0; i < nReps; i++) {
        reroll_HNK(&hnk, bst, 0); // never do rank, it isn't used
        uint32_t gopt = HNK_find_fully_expressed_gopt(&hnk);
        float pj = HNK_calc_p_1(&hnk, backtracking);
        if (HNK_check_gopt(&hnk, gopt)) {
            (*fexn)++;
            (*p_1F) += pj;
        }
        else {
            (*nexn)++;
            (*p_1N) += pj;
        }
    }

    // correct denominators for averages, note that if fexn or nexn is 0 this will return NaN
    (*p_1F) /= (*fexn);
    (*p_1N) /= (*nexn);
    free(backtracking);
}

// determines the average number of local optima and average probability of 1-accessibility
// over nReps NK landscapes with N, K, and the chosen B and K styles
void test_mod_NK(uint32_t N, uint32_t K, B_styles bst, K_styles kst, uint32_t nReps, uint32_t* nopt, float* p_1) {
    NK_landscape nk;
    init_NK(&nk, N, K);
    uint8_t* backtracking = (uint8_t*)malloc((1<<N)*sizeof(uint8_t));
    (*nopt) = 0;
    (*p_1) = 0.0;

    for (uint32_t i = 0; i < nReps; i++) {
        reroll_NK(&nk, kst, bst, 1); // always do rank, several helper functions do not work without it
        (*nopt) += NK_count_local_optima(&nk);
        (*p_1) += NK_calc_p_1(&nk, backtracking);
    }

    (*nopt) /= nReps;
    (*p_1) /= nReps;
    free(backtracking);
}

// perform all of the main experimental cases and print the outputs
void run_tests() {
    //////////////////////////
    // IDEA 1: modified NK
    //////////////////////////
    uint32_t N = 6;
    uint32_t K = 3;

    //////////////////////////
    // Adjacent neighborhoods
    B_styles bst = ADJACENT;

    // Classic NK
    K_styles kst = CLASSIC;

    uint32_t nReps = 10000;
    uint32_t nopt;
    float p_1;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: ADJACENT, CLASSIC; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 1
    kst = RMF_1;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: ADJACENT, RMF theta=1; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 2
    kst = RMF_2;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: ADJACENT, RMF theta=2; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 3
    kst = RMF_3;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: ADJACENT, RMF theta=3; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    //////////////////////////
    // Blocked neighborhoods
    bst = BLOCKED;

    // Classic NK
    kst = CLASSIC;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: BLOCKED, CLASSIC; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 1
    kst = RMF_1;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: BLOCKED, RMF theta=1; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 2
    kst = RMF_2;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: BLOCKED, RMF theta=2; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 3
    kst = RMF_3;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: BLOCKED, RMF theta=3; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    //////////////////////////
    // Random neighborhoods
    bst = RANDOM;
    nReps = 1000000;

    // Classic NK
    kst = CLASSIC;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: RANDOM, CLASSIC; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 1
    kst = RMF_1;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: RANDOM, RMF theta=1; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 2
    kst = RMF_2;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: RANDOM, RMF theta=2; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 3
    kst = RMF_3;
    test_mod_NK(N, K, bst, kst, nReps, &nopt, &p_1);
    printf("IDEA 1: RANDOM, RMF theta=3; nReps = %d; nopt = %d, p_1 = %g\n", nReps, nopt, p_1);

    //////////////////////////
    // IDEA 2: HNK
    //////////////////////////
    N = 8;
    K = 4;
    nReps = 10000;
    uint32_t fexn;
    uint32_t nexn;
    float p_1F;
    float p_1N;

    //////////////////////////
    // H = 1
    uint32_t H = 1;

    // Adjacent neighborhoods
    bst = ADJACENT;
    test_HNK(H, N, K, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=1, ADJACENT; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Blocked neighborhoods
    bst = BLOCKED;
    test_HNK(H, N, K, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=1, BLOCKED; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Random neighborhoods
    bst = RANDOM;
    nReps = 10000000;
    test_HNK(H, N, K, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=1, RANDOM; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);
    
    //////////////////////////
    // H = 2
    H = 2;
    nReps = 10000;

    // Adjacent neighborhoods
    bst = ADJACENT;
    test_HNK(H, N, K, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=2, ADJACENT; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Blocked neighborhoods
    bst = BLOCKED;
    test_HNK(H, N, K, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=2, BLOCKED; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Random neighborhoods
    bst = RANDOM;
    nReps = 10000000;
    test_HNK(H, N, K, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=2, RANDOM; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    //////////////////////////
    // H = 3
    H = 3;
    nReps = 10000;

    // Adjacent neighborhoods
    bst = ADJACENT;
    test_HNK(H, N, K, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=3, ADJACENT; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Blocked neighborhoods
    bst = BLOCKED;
    test_HNK(H, N, K, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=3, BLOCKED; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Random neighborhoods
    bst = RANDOM;
    nReps = 10000000;
    test_HNK(H, N, K, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=3, RANDOM; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);
}

void pretty_print_node(node* n, NK_landscape* nk, uint8_t bitstring_genotype, uint8_t do_rank) {
    uint32_t s = n->genotype;
    uint32_t N = nk->N;
    uint32_t K = nk->K;
    K_neighborhood* nbs = nk->nbs;
    printf("NODE ");
    if (bitstring_genotype) print_genotype_bitstring(s, N);
    else                    printf("%d", s);
    if (do_rank) {
        printf(" %g %d", n->fitness, n->rank);
    }
    else {
        printf(" %g", n->fitness);
    }
    
    // manual recalculation of fitness
    float recalc_fitness = 0.0;
    for (uint32_t i = 0; i < N; i++) {
        // isolate the relevant loci
        uint32_t subgenotype = 0;
        for (uint32_t ell = 0; ell < K; ell++) {
            subgenotype |= ((s >> (nbs[i].B[ell])) & 1) << ell;
        }
        recalc_fitness += nbs[i].fitnesses[subgenotype];
    }
    
    if (recalc_fitness != n->fitness) {
      printf("***** Fitness calculation error! *****");
    }
    
    printf("\n");
}

void pretty_print_HNK_node(node* n, HNK_landscape* hnk, uint8_t bitstring_genotype, uint8_t do_rank) {
    printf("HNKNODE ");
    uint32_t H = hnk->H;
    uint32_t N = hnk->N;
    uint32_t K = hnk->K;
    NK_landscape* gamma_landscape = hnk->gamma_landscape;
    uint32_t* G = hnk->G;
    uint32_t s = n->genotype;
    uint32_t gamma = HNK_sigma_to_gamma(s, H, G);
    if (bitstring_genotype) {
        print_genotype_bitstring(s, N+H);
        printf(" ");
        print_genotype_bitstring(gamma, N);
    }
    else {
        printf("%d", s);
        printf(" %d", gamma);
    }
    
    if (do_rank) {
        printf(" %g %d %g %d", n->fitness, n->rank, gamma_landscape->nodes[gamma].fitness, gamma_landscape->nodes[gamma].rank);
    }
    else {
        printf(" %g %g", n->fitness, gamma_landscape->nodes[gamma].fitness);
    }
    
    // correctness checks
    if (n->fitness != gamma_landscape->nodes[gamma].fitness) {
        printf("***** Fitness mismatch error! *****");
    }
    
    // manual recalculation of gamma
    uint32_t recalc_gamma = s >> H;
    for (uint32_t i = 0; i < N; i++) {
        for (uint32_t h = 0; h < H; h++) {
            if ((G[h] >> i) & 1) { // if this bit is supposed to interact
                if (!((s >> h) & 1)) { // the control gene is off
                    recalc_gamma &= ~(1 << i); // shut off this position
                    break; // no point looking at the others
                }
            }
        }
    }
    if (gamma != recalc_gamma) {
      printf("***** Gamma masking error! *****");
    }
    
    printf("\n");
}

void pretty_print_bitstring_list(uint32_t* list, uint32_t list_len, uint32_t str_len) {
    char ret[32*list_len];
    char* rp = ret;
    for (uint32_t i = 0; i < list_len-1; i++) {
        print_genotype_bitstring(list[i], str_len);
        printf(" ");
    }
    
    print_genotype_bitstring(list[list_len-1], str_len);
}

void pretty_print_integer_list(uint32_t* list, uint32_t list_len) {
    for (uint32_t i = 0; i < list_len-1; i++) {
        printf("%d ", list[i]);
    }
    printf("%d", list[list_len-1]);
}

void pretty_print_K_neighborhood(K_neighborhood* kn, uint32_t N, uint32_t K) {
    printf("KNEIGHBORHOOD | ");
    pretty_print_integer_list(kn->B, K);
    printf(" |");
    for (uint32_t i = 0; i < (1<<K)-1; i++) {
        printf(" %g ", kn->fitnesses[i]);
    }
    printf(" %g\n",  kn->fitnesses[(1<<K)-1]);
}

void pretty_print_NK(NK_landscape* nk, uint8_t genotype_bitstring) {
    printf("NK %d %d %s %s\n", nk->N, nk->K, kst_to_str(nk->kst), bst_to_str(nk->bst));
    for (uint32_t i = 0; i < (1<<(nk->N)); i++) {
        pretty_print_node(&(nk->nodes[i]), nk, genotype_bitstring, 1);
    }
    for (uint32_t i = 0; i < nk->N; i++) {
        pretty_print_K_neighborhood(&(nk->nbs[i]), nk->N, nk->K);
    }
    printf("!NK\n");
}

void pretty_print_HNK(HNK_landscape* hnk, uint8_t genotype_bitstring) {
    printf("HNK %d %d %d %s ", hnk->H, hnk->N, hnk->K, bst_to_str(hnk->bst));
    pretty_print_integer_list(hnk->G, hnk->H);
    printf("\n"); 
    for (uint32_t i = 0; i < (1<<(hnk->N + hnk->H)); i++) {
        pretty_print_HNK_node(&(hnk->nodes[i]), hnk, genotype_bitstring, 1);
    }
    pretty_print_NK(hnk->gamma_landscape, genotype_bitstring);
    printf("!HNK\n");
}

void HNK_audit_localmax(HNK_landscape* hnk) {
    uint8_t* scratch = (uint8_t*)malloc((1<<(hnk->H+hnk->N))*sizeof(uint8_t));
    uint32_t total_maxima = HNK_count_local_optima(hnk, scratch);
    printf("Auditing local maxima (expecting %d):\n", total_maxima);
    uint32_t count = 0;
    for (uint32_t s = 0; s < (1 << (hnk->N)); s++) {
        uint32_t hs = (s << (hnk->H)) + ((1 << (hnk->H))-1);
        uint8_t is_local_max = 1;
        float sfit = hnk->nodes[hs].fitness;
        pretty_print_HNK_node(&(hnk->nodes[hs]), hnk, 1, 1);
        for (uint32_t m = (1<<(hnk->H)); m < (1 << (hnk->N + hnk->H)); m<<=1) {
            uint32_t mod = hs ^ m;
            printf("... ");
            print_genotype_bitstring(mod, hnk->N + hnk->H);
            float mfit = hnk->nodes[mod].fitness;
            if (mfit > sfit) is_local_max = 0;
            printf(": %g\n", mfit);
        }
        if (is_local_max) {
            printf("... ** local maximum **\n");
            count++;
        }
        printf("\n");
    }
    if (total_maxima != count) {
      printf("***** mismatch optima count: %d expected versus %d counted\n", total_maxima, count);
    }
    else {
      printf("Optima count matched.\n");
    }
}

// very slow alternate method of finding the p_1 value
float HNK_alt_p_1(HNK_landscape* hnk, uint8_t* backtracking) {
    uint32_t count = 0;
    for (uint32_t s = 0; s < (1<<(hnk->N)); s++) {
        uint32_t hs = (s << (hnk->H)) + (1<<(hnk->H))-1;
        if (HNK_has_access(hnk, hs, backtracking)) {
            count++;
            printf("[ALT] genotype "); print_genotype_bitstring(hs, hnk->H+hnk->N); printf(" has access\n");
        }
    }
    
    return (count + 0.0)/(1<<(hnk->N));
}

void HNK_audit_p_1(HNK_landscape* hnk) {
    uint8_t* scratch = (uint8_t*)malloc((1<<(hnk->N+hnk->H))*sizeof(uint8_t));
    clock_t start = clock();
    float standard_p_1 = HNK_calc_p_1(hnk, scratch);
    clock_t stop = clock();
    printf("Result from regular (fast) method: %g (%lf seconds)\n", standard_p_1, (float)(stop - start) / CLOCKS_PER_SEC);
    // fflush(stdout);
    start = clock();
    float alt_p_1 = HNK_alt_p_1(hnk, scratch);
    stop = clock();
    printf("Result from hammd (slower) method: %g (%lf seconds)\n", alt_p_1, (float)(stop - start) / CLOCKS_PER_SEC);
    // fflush(stdout);
    if (standard_p_1 != alt_p_1) {
        printf("***** p_1 calculation mismatch *****\n");
    }
    else {
      printf("p_1 calculation matched.\n");
    }
    free(scratch);
}

int main(int argc, char** argv) {
    // HNK landscape example
    srand(0);
    HNK_landscape hnk;
    uint32_t H = 3;
    uint32_t N = 8;
    uint32_t K = 4;
    init_HNK(&hnk, 3, 8, 4);
    
    for (uint32_t i = 0; i < 1000; i++) {
        reroll_HNK(&hnk, RANDOM, 1); // do rank just so we can test it

        pretty_print_HNK(&hnk, 1);
        fflush(stdout);
        HNK_audit_localmax(&hnk);
        fflush(stdout);
        HNK_audit_p_1(&hnk);
        fflush(stdout);
    }

    return EXIT_SUCCESS;
}