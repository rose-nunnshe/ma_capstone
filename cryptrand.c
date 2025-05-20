#include <cryptrand.h>

uint64_t generate_random_64bit_from_timestamp() {
    //static uint64_t last = 0;
    static uint64_t ret = 0;
    ret = rand();
    //uint64_t starting = rand() ^ last;
    //std::string s;
    //s.append(std::to_string(starting));
    //std::hash<std::string> h;
    //uint64_t ret = h(s);
    //last = ret;
    //printf("returning seed %ld\n", ret);
    return ret;
}