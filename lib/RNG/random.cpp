#include "random.hpp"
#include <iostream>

pcg32_random_t::pcg32_random_t()
{
    std::cout << "random setting start" << std::endl;
    //initial seed setting
    std::random_device dev_rnd;
    uint64_t initstate = dev_rnd();
    uint64_t initseq = dev_rnd();

    //initial rng ajustment
    state = 0U;
    inc = (initseq << 1u) | 1u;
    rnd_integer();
    state += initstate;
    rnd_integer();
    std::cout << "random setting fisish" << std::endl;
}

uint32_t pcg32_random_t::rnd_integer()
{
    uint64_t oldstate = state;
    // Advance internal state
    state = oldstate * 6364136223846793005ULL + (inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

float pcg32_random_t::rnd()
{
    return std::ldexp(rnd_integer(), -32);
}