#include <fitness_landscapes.h>

static uint32_t failure_count = 0;
static uint32_t NK_rebuild_count = 0;

static uint8_t global_verbose = 0;

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
void init_NK(NK_landscape* nk, uint32_t N, uint32_t K, B_styles bst) {
    if (global_verbose) {
      printf("<init_NK>\n"); fflush(stdout);
    }
    // basic assignments
    nk->N = N;
    nk->K = K;
    nk->nodes = (node*)malloc((1 << N)*sizeof(node));
    nk->nbs = (K_neighborhood*)malloc(N*sizeof(K_neighborhood));

    // make neighborhoods
    for (uint32_t n = 0; n < N; n++) {
        K_neighborhood* nb = nk->nbs + n;
        nb->K = K;
        nb->B = (uint32_t*)malloc((nb->K)*sizeof(uint32_t));
        nb->fitnesses = (double*)malloc((1<<(nb->K))*sizeof(double));
    }
    
    if (global_verbose) {
      printf("</init_NK>\n"); fflush(stdout);
    }
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
    double diff1 = ((node*)arg1)->fitness - ((node*)arg2)->fitness;
    if (diff1 < 0) return -1;
    if (diff1 > 0) return 1;
    return 0;
}
int compare2(const void* arg1, const void* arg2) {
    return ((node*)arg1)->genotype - ((node*)arg2)->genotype;
}

double NK_calc_fitness(NK_landscape* nk, uint32_t s, uint8_t verbose) {
    uint32_t N = nk->N;
    double fit = 0.0;
    if (verbose) { 
        printf("Calculating fitness of ");
        print_genotype_bitstring(s, N);
        printf("\n");
    }
    for (uint32_t ell = 0; ell < N; ell++) {
        uint32_t subgenotype = 0;
        uint32_t K = nk->nbs[ell].K; // check this neighborhood's K
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

uint8_t NK_rebuild(NK_landscape* nk, uint8_t do_rank) {
  NK_rebuild_count++;
  uint32_t N = nk->N;
  uint8_t succeeded = 1;
  
  // iterate through the genotypes building their fitness
  for (uint32_t s = 0; s < (1<<N); s++) {
      node* current = (nk->nodes) + s;
      current->genotype = s;
      current->fitness = NK_calc_fitness(nk, s, global_verbose); // don't want verbose
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
              if (global_verbose > 10) {
                  printf("***** Holy God, it's a tie! (%f, genotype ", nk->nodes[s].fitness);
                  print_genotype_bitstring(nk->nodes[s].genotype, nk->N);
                  printf(" with %f, genotype ", nk->nodes[s-1].fitness); 
                  print_genotype_bitstring(nk->nodes[s-1].genotype, nk->N); 
                  printf(") *****\n"); 
              }
              
              // need to interrogate exactly how this happened
              NK_calc_fitness(nk, nk->nodes[s].genotype, 1);
              NK_calc_fitness(nk, nk->nodes[s-1].genotype, 1);
              // if we didn't get an injective function, we just have to try again.
              succeeded = 0; 
              failure_count++;
              break;
          }
      }

      // then, sort the nodes by genotype so they're back in index order
      qsort(nk->nodes, (1<<(nk->N)), sizeof(node), compare2);
  }
  return succeeded;
}

// builds neighborhoods, random fitness values, etc.
// DOES NOT allocate anything (N and K already need to be correctly set up)
void reroll_NK(NK_landscape* nk, double theta, B_styles bst, uint8_t do_rank, uint8_t suppress_landscape_update = 0) {
    if (global_verbose) {
      printf("<reroll_NK>\n"); fflush(stdout);
    }
    uint32_t N = nk->N;
    uint32_t K = nk->K;
    nk->theta = theta;
    nk->bst = bst;
    
    if (bst == BLOCKED) { // precalculate the blocks so that we can do it right, here, once.
      uint32_t cl = 0;
      while ((cl + K) < N) {
        // create a block of size K starting at cl
        for (uint32_t i = cl; i < cl+K; i++) {
          K_neighborhood* nb = nk->nbs + i;
          for (uint32_t k = 0; k < K; k++) {
            nb->B[k] = cl+k;
          }
        }
        
        cl += K;
      }
      
      // residual block (if there is one): just the last remaining loci
      uint32_t residual_size = N-cl;
      for (uint32_t i = cl; i < N; i++) {
        K_neighborhood* nb = nk->nbs + i;
        nb->K = residual_size;
        for (uint32_t k = 0; k < residual_size; k++) {
          nb->B[k] = cl+k;
        }
      }
    }
    
    uint8_t succeeded = 0;
    while (!succeeded) {
        succeeded = 1;
        // build neighborhoods for non-blocked styles
        for (uint32_t n = 0; n < N; n++) {
            K_neighborhood* nb = nk->nbs + n;
            uint32_t K = nb->K; // remember, this might not always be the "goal" K of the overall landscape

            for (uint32_t i = 0; i < K; i++) {
                switch (bst) {
                    case BLOCKED:
                        // don't do anything -- the correct indices were handled above
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
                            uint32_t candidate = genrand_int32() % N;
                            uint32_t j = 0;
                            while (j < i) {
                                if (candidate != (nb->B)[j]) { // no collision at this index
                                    j++;
                                }
                                else {
                                    // collision, have to start over
                                    candidate = genrand_int32() % N;
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

            if (!suppress_landscape_update) {            
              // need to choose a peak genotype
              uint32_t delta = genrand_int32() % (1<<K);

              for (uint32_t s = 0; s < (1 << K); s++) {
                  (nb->fitnesses)[s] = -(nk->theta)*hammd(s, delta) + boxmuller();
              }
            }
        }

        succeeded = NK_rebuild(nk, do_rank);
    }
    // so now, all nodes have a correct fitness and correct rank
    if (global_verbose) {
      printf("</reroll_NK>\n"); fflush(stdout);
    }
}

// Allocates nodes and space for underlying NK landscape
// Note: intolerant of H=0 -- if you want to simulate that, use the NK model 
//       with K_style = CLASSIC
void init_HNK(HNK_landscape* hnk, uint32_t H, uint32_t N, uint32_t K, B_styles bst) {
    if (global_verbose) {
      printf("<init_HNK>\n"); fflush(stdout);
    }
    hnk->H = H;
    hnk->N = N;
    hnk->K = K;

    hnk->G = (uint32_t*)malloc(H*sizeof(uint32_t));
    hnk->nodes = (node*)malloc((1<<(H+N))*sizeof(node));
    hnk->gamma_landscape = (NK_landscape*)malloc(sizeof(NK_landscape));
    
    init_NK(hnk->gamma_landscape, N, K, bst);
    if (global_verbose) {
      printf("</init_HNK>\n"); fflush(stdout);
    }
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
void reroll_HNK(HNK_landscape* hnk, double theta, B_styles bst, uint8_t do_rank) {
    if (global_verbose) {
      printf("<reroll_HNK>\n"); fflush(stdout);
    }
    uint32_t H = hnk->H;
    uint32_t N = hnk->N;
    uint32_t K = hnk->K;
    hnk->bst = bst;

    // create the control relationships
    for (uint32_t h = 0; h < H; h++) {
        uint32_t G_h = 0;
        for (uint32_t i = 0; i < K; i++) {
            uint32_t candidate = 1 << (genrand_int32() % N);
            while (G_h & candidate) { // this bit was picked before
                candidate = 1 << (genrand_int32() % N);
            }
            G_h |= candidate;
        }
        (hnk->G)[h] = G_h;
    }

    // build the underlying classic NK
    reroll_NK(hnk->gamma_landscape, theta, bst, do_rank);

    // figure out what correspondence to give based on masking
    for (uint32_t s = 0; s < (1<<(H+N)); s++) {
        (hnk->nodes)[s].genotype = s;

        // construct gamma
        uint32_t gamma = HNK_sigma_to_gamma(s, H, hnk->G);

        // give this genotype the same fitness and rank as its gamma correspondence
        (hnk->nodes)[s].fitness = hnk->gamma_landscape->nodes[gamma].fitness;
        (hnk->nodes)[s].rank = hnk->gamma_landscape->nodes[gamma].rank; // if do_rank = 0 these are just garbage
    }
    if (global_verbose) {
      printf("</reroll_HNK>\n"); fflush(stdout);
    }
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

uint32_t NK_find_gopt_rankless(NK_landscape* nk) {
  double best_fitness = nk->nodes[0].fitness;
  uint32_t best_s = 0;
  
  for (uint32_t s = 1; s < (1<<(nk->N)); s++) {
    if (nk->nodes[s].fitness > best_fitness) {
      if (global_verbose) {
        printf("s = "), print_genotype_bitstring(s, nk->N); printf(" has better fitness (%1.9f) than previous (s = ", nk->nodes[s].fitness); print_genotype_bitstring(best_s, nk->N); printf(", fitness %1.9f)\n", best_fitness);
      }
      best_fitness = nk->nodes[s].fitness;
      best_s = s;
    }
    else {
      if (global_verbose) {
        printf("s = "), print_genotype_bitstring(s, nk->N); printf(" has worse fitness (%1.9f) than previous (s = ", nk->nodes[s].fitness); print_genotype_bitstring(best_s, nk->N); printf(", fitness %1.9f)\n", best_fitness);
      }
    }
  }
  
  return best_s;
}

uint32_t NK_calc_p_1_recur(NK_landscape* nk, uint8_t* backtracking, uint32_t s, uint32_t depth) {
    backtracking[s] = 1;

    uint32_t ret = 1; // this one
    uint32_t N = nk->N;
    
    if (global_verbose) {
      for (int i = 0; i < depth; i++) {
        for (int j = 0; j < N+16; j++) {
          printf(" ");
        }
      }
      print_genotype_bitstring(nk->nodes[s].genotype, N); printf(" (%2.9f)*", nk->nodes[s].fitness);
    }

    for (uint32_t m = 1; m < (1<<N); m<<=1) {
        uint32_t next = s ^ m;
        if ((nk->nodes[next].fitness < nk->nodes[s].fitness) && (!backtracking[next])) { // fitness must be decreasing working away from global opt 
            if (global_verbose) printf("\n");
            ret += NK_calc_p_1_recur(nk, backtracking, next, depth+1);
        }
        else {
          if (m == (1<<(N-1))) {
            if (global_verbose) {
              printf("** Refusing to recur: think that next fitness %2.9f is higher than this", nk->nodes[next].fitness);
            }
          }
        }
    }
    return ret;
}

// calculates p_1 by working from the global optimum outward
double NK_calc_p_1(NK_landscape* nk, uint8_t* backtracking) {
    for (uint32_t s = 0; s < (1<<(nk->N)); s++) {
        backtracking[s] = 0;
    }

    return ((double)NK_calc_p_1_recur(nk, backtracking, NK_find_gopt_rankless(nk), 0))/(1<<(nk->N));
}

double NK_alt_p_1(NK_landscape* nk, uint8_t* backtracking) {
    uint32_t count = 0;
    for (uint32_t s = 0; s < (1<<(nk->N)); s++) {
        if (NK_has_access(nk, s, backtracking)) count++;
    }
    
    return ((double)count)/(1<<(nk->N));
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
    double max_fit = hnk->nodes[h_mask].fitness;
    uint32_t max_s = hnk->nodes[h_mask].genotype;

    for (uint32_t s = (1<<H) + h_mask; s < (1<<(N+H)); s += (1<<H)) {
        double this_fit = hnk->nodes[s].fitness;
        if (this_fit > max_fit) {
            max_fit = this_fit;
            max_s = s;
        }
    }

    return max_s;
}

uint32_t HNK_calc_p_1_recur(HNK_landscape* hnk, uint8_t* backtracking, uint32_t s, uint32_t recursion_depth) {
    if (backtracking[s]) return 0;
    backtracking[s] = 1;
    if (global_verbose) {
      printf("[HNK_calc_p_1_recur] recursion depth = %d... ", recursion_depth); fflush(stdout);
    }
    uint32_t ret = 0;
    uint32_t H = hnk->H;
    uint32_t N = hnk->N;
    uint32_t h_mask = (1<<H)-1;
    if ((s & h_mask) == h_mask) {
        ret += 1; // this is a fully-expressed genotype
        if (global_verbose) {
          printf("[REGULAR] genotype "); print_genotype_bitstring(s, hnk->H+hnk->N); printf(" has access"); fflush(stdout);
        }
    }
    
    if (global_verbose) {
      printf("\n"); fflush(stdout);
    }
    
    // look for neighbors
    for (uint32_t m = 1; m < (1<<(H+N)); m<<=1) {
        uint32_t next = s ^ m;
        if (hnk->nodes[next].fitness <= hnk->nodes[s].fitness) { // fitness must be non-increasing when moving away from optimum
            ret += HNK_calc_p_1_recur(hnk, backtracking, next, recursion_depth+1);
        }
    }

    return ret;
}

double HNK_calc_p_1_norecur(HNK_landscape* hnk, uint8_t* backtracking, uint32_t gopt) {
  // instead of recurrence, we maintain our own "stack" dedicated to traversing the landscape
  // this function is what forced us to go to C++ so that I could use STL vector rather than write my own
  double ret = 0;
  std::vector<uint32_t> callstack = std::vector<uint32_t>();
  std::vector<uint32_t> mstack = std::vector<uint32_t>();
  callstack.push_back(gopt);
  mstack.push_back(1);
  uint32_t callstack_size = 1;
  
  uint32_t H = hnk->H;
  uint32_t N = hnk->N;
  uint32_t h_mask = (1<<H)-1;
  
  while (callstack_size > 0) {
    // get the callstack head
    uint32_t s = callstack[callstack_size-1];
    uint32_t m = mstack[callstack_size-1];
    if (global_verbose) {
      printf("[HNK_calc_p_1_norecur] s = "); print_genotype_bitstring(s, H+N); printf(" (%d) m = %d callstack_size = %d backtracking = %d\n", s, m, callstack_size, backtracking[s]);
    }
    
    if (!backtracking[s]) { // only on the first time through, allow for this genotype to be counted
      if ((s & h_mask) == h_mask) {
        ret += 1; // this is a fully-expressed genotype
        if (global_verbose) {
          printf("[HNK_calc_p_1_norecur] genotype "); print_genotype_bitstring(s, hnk->H+hnk->N); printf(" has access\n"); fflush(stdout);
        }
      }
    }
    
    backtracking[s] = 1;
    
    // look for neighbors
    if (m < (1<<(H+N))) {
      uint32_t next = s ^ m;
      mstack[callstack_size-1] <<= 1;
      if (hnk->nodes[next].fitness <= hnk->nodes[s].fitness) { // fitness must be non-increasing when moving away from optimum
        if (!backtracking[next]) {
          callstack.push_back(next);
          mstack.push_back(1); // for the next node
          callstack_size++;
        }
      }
    }
    else {
      // this node is done
      callstack.pop_back();
      mstack.pop_back();
      callstack_size--;
    }
  }
  
  return ret/(1<<N);
}

// calculate p_1 for a given HNK landscape (requires rank data)
double HNK_calc_p_1(HNK_landscape* hnk, uint8_t* backtracking) {
    // set up anti-backtracking helpers
    for (uint32_t i = 0; i < (1<<(hnk->N+hnk->H)); i++) {
        backtracking[i] = 0;
    }

    uint32_t gopt = HNK_find_fully_expressed_gopt(hnk);
    #ifdef USE_RECURSION
      // elaborate from gopt with non-increasing fitness, only counting fully expressed genotypes
      // but being allowed to visit non-fully expressed genotypes during recursion
      double ret = ((double)HNK_calc_p_1_recur(hnk, backtracking, gopt, 0))/(1<<(hnk->N));
    #else
      double ret = HNK_calc_p_1_norecur(hnk, backtracking, gopt);
    #endif
    if (global_verbose) {
      printf("[HNK_calc_p_1] finished returning %f\n", ret); fflush(stdout);
    }
    return ret;
}

void pretty_print_node(node* n, NK_landscape* nk, uint8_t bitstring_genotype, uint8_t do_rank) {
    uint32_t s = n->genotype;
    uint32_t N = nk->N;
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
    double recalc_fitness = 0.0;
    for (uint32_t i = 0; i < N; i++) {
        // isolate the relevant loci
        uint32_t subgenotype = 0;
        uint32_t K = nbs[i].K;
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

void pretty_print_K_neighborhood(K_neighborhood* kn, uint32_t N) {
    uint32_t K = kn->K;
    printf("KNEIGHBORHOOD | ");
    pretty_print_integer_list(kn->B, K);
    printf(" |");
    for (uint32_t i = 0; i < (1<<K)-1; i++) {
        printf(" %g ", kn->fitnesses[i]);
    }
    printf(" %g\n",  kn->fitnesses[(1<<K)-1]);
}

void pretty_print_NK(NK_landscape* nk, uint8_t genotype_bitstring) {
    printf("NK %d %d %f %s\n", nk->N, nk->K, (nk->theta), bst_to_str(nk->bst));
    for (uint32_t i = 0; i < (1<<(nk->N)); i++) {
        pretty_print_node(&(nk->nodes[i]), nk, genotype_bitstring, 1);
    }
    for (uint32_t i = 0; i < nk->N; i++) {
        pretty_print_K_neighborhood(&(nk->nbs[i]), nk->N);
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

uint32_t HNK_audit_localmax(HNK_landscape* hnk) {
    uint8_t* scratch = (uint8_t*)malloc((1<<(hnk->H+hnk->N))*sizeof(uint8_t));
    uint32_t total_maxima = HNK_count_local_optima(hnk, scratch);
    printf("Auditing local maxima (expecting %d):\n", total_maxima);
    uint32_t count = 0;
    for (uint32_t s = 0; s < (1 << (hnk->N)); s++) {
        uint32_t hs = (s << (hnk->H)) + ((1 << (hnk->H))-1);
        uint8_t is_local_max = 1;
        double sfit = hnk->nodes[hs].fitness;
        pretty_print_HNK_node(&(hnk->nodes[hs]), hnk, 1, 1);
        for (uint32_t m = (1<<(hnk->H)); m < (1 << (hnk->N + hnk->H)); m<<=1) {
            uint32_t mod = hs ^ m;
            printf("... ");
            print_genotype_bitstring(mod, hnk->N + hnk->H);
            double mfit = hnk->nodes[mod].fitness;
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
    
    return count;
}

uint32_t NK_audit_localmax(NK_landscape* nk) {
    uint32_t total_maxima = NK_count_local_optima(nk);
    printf("Auditing local maxima (expecting %d):\n", total_maxima);
    uint32_t count = 0;
    for (uint32_t s = 0; s < (1 << (nk->N)); s++) {
        uint8_t is_local_max = 1;
        double sfit = nk->nodes[s].fitness;
        pretty_print_node(&(nk->nodes[s]), nk, 1, 1);
        for (uint32_t m = 1; m < (1 << (nk->N)); m<<=1) {
            uint32_t mod = s ^ m;
            printf("... ");
            print_genotype_bitstring(mod, nk->N);
            double mfit = nk->nodes[mod].fitness;
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
    
    return count;
}

// very slow alternate method of finding the p_1 value
double HNK_alt_p_1(HNK_landscape* hnk, uint8_t* backtracking) {
    uint32_t count = 0;
    for (uint32_t s = 0; s < (1<<(hnk->N)); s++) {
        uint32_t hs = (s << (hnk->H)) + (1<<(hnk->H))-1;
        if (HNK_has_access(hnk, hs, backtracking)) {
            count++;
            // printf("[ALT] genotype "); print_genotype_bitstring(hs, hnk->H+hnk->N); printf(" has access\n");
        }
    }
    
    return ((double)count)/(1<<(hnk->N));
}

double HNK_audit_p_1(HNK_landscape* hnk) {
    uint8_t* scratch = (uint8_t*)malloc((1<<(hnk->N+hnk->H))*sizeof(uint8_t));
    clock_t start = clock();
    double standard_p_1 = HNK_calc_p_1(hnk, scratch);
    clock_t stop = clock();
    printf("Result from regular (fast) method: %g (%lf seconds)\n", standard_p_1, (double)(stop - start) / CLOCKS_PER_SEC);
    // fflush(stdout);
    start = clock();
    double alt_p_1 = HNK_alt_p_1(hnk, scratch);
    stop = clock();
    printf("Result from counting (slower) method: %g (%lf seconds)\n", alt_p_1, (double)(stop - start) / CLOCKS_PER_SEC);
    // fflush(stdout);
    if (standard_p_1 != alt_p_1) {
        printf("***** p_1 calculation mismatch *****\n");
    }
    else {
      printf("p_1 calculation matched.\n");
    }
    free(scratch);
    
    return standard_p_1;
}

double NK_audit_p_1(NK_landscape* nk) {
    uint8_t* scratch = (uint8_t*)malloc((1<<(nk->N))*sizeof(uint8_t));
    clock_t start = clock();
    double standard_p_1 = NK_calc_p_1(nk, scratch);
    clock_t stop = clock();
    printf("Result from regular (fast) method: %g (%lf seconds)\n", standard_p_1, (double)(stop - start) / CLOCKS_PER_SEC);
    // fflush(stdout);
    start = clock();
    double alt_p_1 = NK_alt_p_1(nk, scratch);
    stop = clock();
    printf("Result from counting (slower) method: %g (%lf seconds)\n", alt_p_1, (double)(stop - start) / CLOCKS_PER_SEC);
    // fflush(stdout);
    if (standard_p_1 != alt_p_1) {
        printf("***** p_1 calculation mismatch *****\n");
    }
    else {
      printf("p_1 calculation matched.\n");
    }
    free(scratch);
    
    return standard_p_1;
}

// will require optima to be allocated larger than the actual number turns out
void NK_list_optima(NK_landscape* nk, uint32_t* optima, uint32_t* nOpt) {
  uint32_t count = 0;
  
  for (uint32_t s = 0; s < (1<<(nk->N)); s++) {
    uint8_t opt = 1;
    for (uint32_t m = 1; m < (1<<(nk->N)); m<<=1) {
      if (nk->nodes[s ^ m].fitness > nk->nodes[s].fitness) {
        opt = 0;
        break;
      }
    }
    if (opt) {
      optima[count++] = s;
    }
  }
  
  *nOpt = count;
}

double NK_calc_normalized_rank(NK_landscape* nk, uint32_t s) {
  uint32_t gt = 0;
  const double num_neighbors = nk->N;
  
  for (uint32_t m = 1; m < (1<<(nk->N)); m<<=1) {
    if (nk->nodes[s ^ m].fitness < nk->nodes[s].fitness) {
      gt++;
    }
  }
  
  return ((double)gt)/num_neighbors;
}

int compare3(const void* arg1, const void* arg2) {
    double diff = (((double*)arg1) - ((double*)arg2));
    if (diff < 0) return -1;
    if (diff > 0) return 1;
    else          return 0;
}

void NK_calc_D(NK_landscape* nk, uint32_t* optima, uint32_t nOpt, double* D) {
  double nranks[nOpt];
  for (uint32_t i = 0; i < nOpt; i++) {
    nranks[i] = NK_calc_normalized_rank(nk, optima[i]);
  }
  
  qsort(nranks, nOpt, sizeof(double), compare3);
  /*
  printf("# sorted nranks, nOpt = %d\n# ", nOpt);
  for (uint32_t i = 0; i < nOpt; i++) {
    printf("%f ", nranks[i]);
  }
  printf("\n");
  */
  
  D[0] = 0.0;
  double prev = 0.0;
  uint32_t nri = 0;
  for (uint32_t i = 1; i < D_PRECISION; i++) {
    D[i] = prev;
    while ((nri < nOpt) && (i >= nranks[nri]*(D_PRECISION-1))) { // has to be while because it's possible to have ties
      D[i] = ((double)(nri+1))/nOpt;
      nri++;
      prev = D[i];
    }
  }
}

void NK_inject_one_disruption(NK_landscape* nk, double intensity) {
  uint32_t neighborhood = genrand_int32() % (nk->N);
  uint32_t nbK = nk->nbs[neighborhood].K;
  uint32_t subgenotype = genrand_int32() % (1<<(nbK));
  
  nk->nbs[neighborhood].fitnesses[subgenotype] += intensity*boxmuller();
  
  NK_rebuild(nk, 1); // update ranks for the sake of not making a mess, it might be possible to skip
}

void NK_inject_skew_disruption(NK_landscape* nk, double intensity) {
  uint32_t target_neighborhood = genrand_int32() % (nk->N);
  uint32_t nbK = nk->nbs[target_neighborhood].K;
  uint32_t target_subgenotype = genrand_int32() % (1<<(nbK));
  
  double skew = intensity*boxmuller();
  
  for (uint32_t interactor = 0; interactor < (nbK); interactor++) {
    uint32_t current_locus = nk->nbs[target_neighborhood].B[interactor];
    for (uint32_t subg = 0; subg < (1<<(nbK)); subg++) {
      double multiplier = -1.0;
      if (((subg >> interactor) & 1) == ((target_subgenotype >> interactor) & 1)) {
        multiplier = 1.0;
      }
      if (current_locus != target_neighborhood) {
        multiplier /= (nbK);
      }
      nk->nbs[current_locus].fitnesses[subg] += multiplier*skew;
    }
  }
  
  NK_rebuild(nk, 1);
}

void NK_inject_vector_disruption(NK_landscape* nk, double intensity) {
  double skew[nk->N];
  for (uint32_t i = 0; i < (nk->N); i++) {
    skew[i] = intensity*boxmuller();
  }
  
  for (uint32_t nb = 0; nb < (nk->N); nb++) {
    K_neighborhood* current_nb = &(nk->nbs[nb]);
    uint32_t K = current_nb->K;
    for (uint32_t sg = 0; sg < (1<<(K)); sg++) {
      for (uint32_t b = 0; b < (K); b++) {
        int32_t locus_value = (sg >> b) & 1;
        current_nb->fitnesses[sg] += (2*locus_value-1)*skew[current_nb->B[b]];
      }
    }
  }
  
  NK_rebuild(nk, 1);
}

void print_D(double* D, uint32_t dis_count) {
  for (uint32_t i = 0; i < D_PRECISION; i++) {
    printf("%d %f %f\n", dis_count, ((double)i)/(D_PRECISION-1), D[i]);
  }
  printf("\n");
}

// over nReps HNK landscapes, count how many had a fully expressed global optimum, how many had a 
// non-fully expressed global optimum, and p_1 for each case
void test_HNK(uint32_t H, uint32_t N, uint32_t K, double theta, B_styles bst, uint32_t nReps, uint32_t* fexn, uint32_t* nexn, double* p_1F, double* p_1N) {
    HNK_landscape hnk;
    init_HNK(&hnk, H, N, K, bst);
    uint8_t* backtracking = (uint8_t*)malloc((1<<(H+N))*sizeof(uint8_t));
    (*fexn) = 0;
    (*nexn) = 0;
    (*p_1F) = 0.0;
    (*p_1N) = 0.0;
    for (uint32_t i = 0; i < nReps; i++) {
        reroll_HNK(&hnk, theta, bst, 0);
        uint32_t gopt = HNK_find_fully_expressed_gopt(&hnk);
        double pj = HNK_calc_p_1(&hnk, backtracking);
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
    if (*fexn == 0) {
      *p_1F = 0;
    }
    if (*nexn == 0) {
      *p_1N = 0;
    }
    free(backtracking);
    deinit_HNK(&hnk);
}

// determines the average number of local optima and average probability of 1-accessibility
// over nReps NK landscapes with N, K, and the chosen B and K styles
void test_mod_NK(uint32_t N, uint32_t K, B_styles bst, double theta, uint32_t nReps, double* nopt, double* p_1, uint8_t do_rank = 0) {
    NK_landscape nk;
    init_NK(&nk, N, K, bst);
    uint8_t* backtracking = (uint8_t*)malloc((1<<N)*sizeof(uint8_t));
    (*nopt) = 0.0;
    (*p_1) = 0.0;

    for (uint32_t i = 0; i < nReps; i++) {
        reroll_NK(&nk, theta, bst, do_rank); // do not need rank
        (*nopt) += NK_count_local_optima(&nk);
        (*p_1) += NK_calc_p_1(&nk, backtracking);
        
        if (global_verbose) {
          pretty_print_NK(&nk, 1);
        }
    }

    (*nopt) /= nReps;
    (*p_1) /= nReps;
    free(backtracking);
    deinit_NK(&nk);
}

// same as test_mod_NK, but attempts to estimate the right number of repetitions to run by looking 
// at a rough estimate of precision
void adaptive_test_mod_NK(uint32_t N, uint32_t K, B_styles bst, double theta, uint32_t min_reps, uint32_t max_reps, double desired_precision, double* nopt, double* p_1, double* nReps) {
  NK_landscape nk;
  init_NK(&nk, N, K, bst);
  uint8_t* backtracking = (uint8_t*)malloc((1<<N)*sizeof(uint8_t));
  (*nopt) = 0.0;
  (*p_1) = 0.0;
  
  uint8_t special_case_flag = 0;
  if (bst == BLOCKED) {
    if (((N-1)/K)*K+1 == N) {
      if (N>1) {
          // K divides N-1, so we should replicate all results with masking the last
          special_case_flag = 1;
      }
    }
  }
  
  for (uint32_t i = 0; i < min_reps; i++) {
	uint32_t seed = generate_random_64bit_from_timestamp() % (1ULL << 31); 
	init_genrand(seed); // should prevent trials from influencing each other
  reroll_NK(&nk, theta, bst, 0); // do not need rank
  
  double this_nopt = NK_count_local_optima(&nk);
  (*nopt) += this_nopt;
  double this_p_1 = NK_calc_p_1(&nk, backtracking);
  (*p_1) += this_p_1;
  
  if (special_case_flag) {
    init_genrand(seed); // restore the state so that we will get the exact same sublandscape
    nk.N = nk.N-1;
    reroll_NK(&nk, theta, bst, 0);
    double masked_nopt = NK_count_local_optima(&nk);
    double masked_p_1 = NK_calc_p_1(&nk, backtracking);
    if ((masked_nopt != this_nopt) || (masked_p_1 != this_p_1)) {
      printf("**** Fatal error has occurred. ****\n"); fflush(stdout);
      global_verbose = 1;
      printf("**** Masked landscape: ****\n");
      init_genrand(seed);
      reroll_NK(&nk, theta, bst, 1);
      pretty_print_NK(&nk, 1);
      NK_count_local_optima(&nk);
      NK_calc_p_1(&nk, backtracking);      
      
      nk.N = nk.N+1;
      printf("**** Unmasked landscape: ****\n");
      init_genrand(seed);
      reroll_NK(&nk, theta, bst, 1);
      pretty_print_NK(&nk, 1);
      NK_count_local_optima(&nk);
      NK_calc_p_1(&nk, backtracking);     
      fflush(stdout);        
      exit(0); // just stop
    }
    nk.N = nk.N+1;
    }
  }
  
  if (global_verbose) {
    pretty_print_NK(&nk, 1);
  }
  
  uint32_t total_reps = min_reps;
  
  // perform additional reps waiting for either the max_reps limit to hit or for 
  // the precision to appear to be arrived at
  uint32_t streak = 0;
  while (total_reps < max_reps) {
    reroll_NK(&nk, theta, bst, 0);
    total_reps++;
    
    double next_optc = NK_count_local_optima(&nk);
    double next_p1 = NK_calc_p_1(&nk, backtracking);
    
    double nopt2 = (*nopt) + next_optc;
    double p_1_2 = (*p_1) + next_p1;
    
    // these are always non-negative because we used the same denominator (not seen here b/c total_reps multiplies out)
    double nopt_prec = 2*(nopt2 - (*nopt))/(nopt2 + (*nopt));
    double p1_prec = 2*(p_1_2 - (*p_1))/(p_1_2 + (*p_1));
    
    (*nopt) = nopt2;
    (*p_1) = p_1_2;
    
    if ((nopt_prec < desired_precision) && (p1_prec < desired_precision)) {
      streak++;
      if (streak > 0.1*min_reps) {
        break;
      }
    }
  }
  
  (*nopt) /= total_reps;
  (*p_1) /= total_reps;
  (*nReps) = total_reps;
  free(backtracking);
  deinit_NK(&nk);
}

// perform all of the main experimental cases and print the outputs
void run_tests() {
    //////////////////////////
    // IDEA 1: modified NK
    //////////////////////////
    uint32_t N = 9;
    uint32_t K = 3;

    //////////////////////////
    // Adjacent neighborhoods
    B_styles bst = ADJACENT;

    // Classic NK
    double theta = 0;

    uint32_t nReps = 1000000;
    double nopt;
    double p_1;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: ADJACENT, CLASSIC; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 1
    theta = 1;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: ADJACENT, RMF theta=1; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 2
    theta = 2;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: ADJACENT, RMF theta=2; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 3
    theta = 3;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: ADJACENT, RMF theta=3; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    //////////////////////////
    // Blocked neighborhoods
    bst = BLOCKED;

    // Classic NK
    theta = 0;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: BLOCKED, CLASSIC; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 1
    theta = 1;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: BLOCKED, RMF theta=1; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 2
    theta = 2;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: BLOCKED, RMF theta=2; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 3
    theta = 3;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: BLOCKED, RMF theta=3; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    //////////////////////////
    // Random neighborhoods
    bst = RANDOM;
    nReps = 100000000;

    // Classic NK
    theta = 0;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: RANDOM, CLASSIC; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 1
    theta = 1;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: RANDOM, RMF theta=1; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 2
    theta = 2;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: RANDOM, RMF theta=2; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    // RMF theta = 3
    theta = 3;
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1);
    printf("IDEA 1: RANDOM, RMF theta=3; nReps = %d; nopt = %f, p_1 = %g\n", nReps, nopt, p_1);

    //////////////////////////
    // IDEA 2: HNK
    //////////////////////////
    N = 9;
    K = 3;
    nReps = 1000000;
    uint32_t fexn;
    uint32_t nexn;
    double p_1F;
    double p_1N;

    //////////////////////////
    // H = 1
    uint32_t H = 1;

    // Adjacent neighborhoods
    bst = ADJACENT;
    test_HNK(H, N, K, 0, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=1, ADJACENT; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Blocked neighborhoods
    bst = BLOCKED;
    test_HNK(H, N, K, 0, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=1, BLOCKED; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Random neighborhoods
    bst = RANDOM;
    nReps = 100000000;
    test_HNK(H, N, K, 0, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=1, RANDOM; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);
    
    //////////////////////////
    // H = 2
    H = 2;
    nReps = 1000000;

    // Adjacent neighborhoods
    bst = ADJACENT;
    test_HNK(H, N, K, 0, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=2, ADJACENT; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Blocked neighborhoods
    bst = BLOCKED;
    test_HNK(H, N, K, 0, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=2, BLOCKED; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Random neighborhoods
    bst = RANDOM;
    nReps = 100000000;
    test_HNK(H, N, K, 0, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=2, RANDOM; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    //////////////////////////
    // H = 3
    H = 3;
    nReps = 1000000;

    // Adjacent neighborhoods
    bst = ADJACENT;
    test_HNK(H, N, K, 0, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=3, ADJACENT; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Blocked neighborhoods
    bst = BLOCKED;
    test_HNK(H, N, K, 0, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=3, BLOCKED; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);

    // Random neighborhoods
    bst = RANDOM;
    nReps = 100000000;
    test_HNK(H, N, K, 0, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
    printf("IDEA 2: H=3, RANDOM; nReps = %d; fexn = %d, nexn = %d, p_1F = %g, p_1N = %g\n", nReps, fexn, nexn, p_1F, p_1N);
}

void run_disrupt_test() {
  NK_landscape nk;
  uint32_t N = 12;
  uint32_t K = 4;
  double theta = 0;
  B_styles bst = RANDOM;
  init_NK(&nk, N, K, bst);
  reroll_NK(&nk, theta, bst, 1);
  
  double* D = (double*)malloc(D_PRECISION*sizeof(double));
  uint32_t* optima = (uint32_t*)malloc(D_PRECISION*sizeof(uint32_t));
  uint32_t nOpt = 0;
  NK_list_optima(&nk, optima, &nOpt);
  NK_calc_D(&nk, optima, nOpt, D);
  print_D(D, 0); fflush(stdout);
  // printf("%d %d\n", 0, NK_count_local_optima(&nk));
  
  uint32_t num_disrupts = 100000;
  double intensity = 0.01;
  uint32_t decimate = 100;
  for (uint32_t i = 0; i < num_disrupts; i++) {
    // NK_inject_one_disruption(&nk, intensity); // old style, one random change in one place
    NK_inject_skew_disruption(&nk, intensity); // new style, skews groups of genotypes
    // NK_inject_vector_disruption(&nk, intensity); // newer style, assigns skews per locus
    
    NK_calc_D(&nk, optima, nOpt, D);
    if ((i % decimate) == 0) {
      print_D(D, i+1); fflush(stdout);
      // printf("%d %d\n", i+1, NK_count_local_optima(&nk));
    }
    // printf("Run %d of %d\n", i, num_disrupts); fflush(stdout);
  }
  
  fflush(stdout);
  deinit_NK(&nk);
  free(D);
  free(optima);
}

void run_sweep_test() {
  B_styles bst = BLOCKED;
  
  uint32_t min_N = 20;
  uint32_t max_N = 20;
  uint32_t min_K = 14;
  
  uint32_t min_reps = 3000;
  uint32_t max_reps = 3000;
  double desired_precision = 0.001;
  double nopt;
  double p_1;
  double nReps;
  
  for (uint32_t N = min_N; N <= max_N; N++) {
    for (uint32_t K = min_K; K <= N; K++) {
    //for (uint32_t K = N; K <= N; K++) {
      double theta = 0;
      adaptive_test_mod_NK(N, K, bst, theta, min_reps, max_reps, desired_precision, &nopt, &p_1, &nReps);
      printf("%d %d %f %f %f %f\n", N, K, theta, nopt, p_1, nReps); fflush(stdout);
      theta = 1.0;
      adaptive_test_mod_NK(N, K, bst, theta, min_reps, max_reps, desired_precision, &nopt, &p_1, &nReps);
      printf("%d %d %f %f %f %f\n", N, K, theta, nopt, p_1, nReps); fflush(stdout);
      theta = 2.0;
      adaptive_test_mod_NK(N, K, bst, theta, min_reps, max_reps, desired_precision, &nopt, &p_1, &nReps);
      printf("%d %d %f %f %f %f\n", N, K, theta, nopt, p_1, nReps); fflush(stdout);
      theta = 3.0;
      adaptive_test_mod_NK(N, K, bst, theta, min_reps, max_reps, desired_precision, &nopt, &p_1, &nReps);
      printf("%d %d %f %f %f %f\n", N, K, theta, nopt, p_1, nReps); fflush(stdout);
      //theta = 10.0;
      //adaptive_test_mod_NK(N, K, bst, theta, min_reps, max_reps, desired_precision, &nopt, &p_1, &nReps);
      //printf("%d %d %f %f %f %f\n", N, K, theta, nopt, p_1, nReps); fflush(stdout);
    }
    printf("\n"); fflush(stdout);
  }
}

void run_HNK_sweep_test() {
  B_styles bst = ADJACENT;
  
  uint32_t min_N = 20;
  uint32_t max_N = 20;
  uint32_t min_K = 6;
  
  uint32_t min_reps = 3000;
  uint32_t max_reps = 3000;
  double desired_precision = 0.001;
  uint32_t nexn;
  uint32_t fexn;
  double p_1F;
  double p_1N;
  uint32_t nReps = 3000;
  
  for (uint32_t N = min_N; N <= max_N; N++) {
    for (uint32_t K = min_K; K <= N; K++) {
      for (double theta = 0; theta <= 0; theta += 1) {
        for (uint32_t H = 2; H < 3; H++) {
          test_HNK(H, N, K, theta, bst, nReps, &fexn, &nexn, &p_1F, &p_1N);
          printf("%d %d %d %f %d %d %f %f\n", H, N, K, theta, fexn, nexn, p_1F, p_1N); fflush(stdout);
        }
      }
    }
    printf("\n"); fflush(stdout);
  }
}

void run_interpblock_test() {
  NK_landscape nk;
  uint32_t N = 15;
  uint32_t K = 5;

  B_styles bst = BLOCKED;
  init_NK(&nk, N, K, bst);
  
  // Classic NK
  double theta = 0;

  reroll_NK(&nk, theta, bst, 1);
  pretty_print_NK(&nk, 1);
  NK_audit_localmax(&nk);
  NK_audit_p_1(&nk);
  
  deinit_NK(&nk);
}

void run_single_test() {
  NK_landscape nk;
  
  B_styles bst = RANDOM;
  
  uint32_t N = 8;
  uint32_t K = 4;
  
  uint32_t min_reps = 1;
  uint32_t max_reps = 1;
  double desired_precision = 0.001;
  double nopt;
  double p_1;
  double nReps;
  
  init_NK(&nk, N, K, bst);
  double theta = 0;
  adaptive_test_mod_NK(N, K, bst, theta, min_reps, max_reps, desired_precision, &nopt, &p_1, &nReps);
  printf("[N = %d][K = %d][CLASSIC][nReps = %f] nopt = %f, p_1 = %f\n", N, K, nReps, nopt, p_1);
  theta = 1;
  adaptive_test_mod_NK(N, K, bst, theta, min_reps, max_reps, desired_precision, &nopt, &p_1, &nReps);
  printf("[N = %d][K = %d][RMF_1][nReps = %f] nopt = %f, p_1 = %f\n", N, K, nReps, nopt, p_1);
  theta = 2;
  adaptive_test_mod_NK(N, K, bst, theta, min_reps, max_reps, desired_precision, &nopt, &p_1, &nReps);
  printf("[N = %d][K = %d][RMF_2][nReps = %f] nopt = %f, p_1 = %f\n", N, K, nReps, nopt, p_1);
  theta = 3;
  adaptive_test_mod_NK(N, K, bst, theta, min_reps, max_reps, desired_precision, &nopt, &p_1, &nReps);
  printf("[N = %d][K = %d][RMF_3][nReps = %f] nopt = %f, p_1 = %f\n", N, K, nReps, nopt, p_1);
  deinit_NK(&nk);
}

void run_mystery_test() {
  // BROKEN, DO NOT RUN
  uint32_t N = 11;
  uint32_t K = 2;
  B_styles bst = BLOCKED;
  double theta = 0.0;
  uint32_t nReps = 1;
  double nopt = 0;
  double p_1 = 0;
  
  for (int i = 1; i < 1000000; i++) {
    //init_genrand(i);
    test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1, 1);
    double p_1_1 = p_1;
    printf("Without mask: nopt = %f p_1 = %f\n", nopt, p_1);
    //init_genrand(i);
    test_mod_NK(N-1, K, bst, theta, nReps, &nopt, &p_1, 1);
    printf("With mask: nopt = %f p_1 = %f\n", nopt, p_1);
    
    if (p_1_1 != p_1) {
      printf("***** Error found! i = %d *****\n", i);
      global_verbose = 1;
      //init_genrand(i);
      test_mod_NK(N, K, bst, theta, nReps, &nopt, &p_1, 1);
      printf("Without mask: nopt = %f p_1 = %f\n", nopt, p_1);
      //init_genrand(i);
      test_mod_NK(N-1, K, bst, theta, nReps, &nopt, &p_1, 1);
      printf("With mask: nopt = %f p_1 = %f\n", nopt, p_1);
      break;
    }
  }
}

int main(int argc, char** argv) {
    //init_genrand(0xa3a5b7d8UL); // this number has no significance, it just has to be a 32-bit, non-zero integer
    failure_count = 0;
    NK_rebuild_count = 0;
    
    clock_t start = clock();
    // run_tests();
    // run_disrupt_test();
    // run_sweep_test();
    // run_mystery_test();
    run_HNK_sweep_test();
    // run_interpblock_test();
    // run_single_test();
    clock_t stop = clock();
    double time_s = (double)(stop - start) / CLOCKS_PER_SEC;
    printf("# Final time: %d h : %d m : %d s\n", (uint64_t)(time_s/3600), (uint64_t)((time_s - 3600*(uint64_t)(time_s/3600))/60), (uint64_t)(time_s - 60*(uint64_t)(time_s/60)));
    
    printf("# Final NK failure count: %d out of a total of %d rerolls (%g%%)\n", failure_count, NK_rebuild_count, (double)failure_count/NK_rebuild_count*100);

    return EXIT_SUCCESS;
}