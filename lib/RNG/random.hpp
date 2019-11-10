#ifndef RANDOM_HPP
#define RANDOM_HPP
#include <cstdint>
#include <cmath>
#include <random>

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
{
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
}

void initRNG(pcg32_random_t* rng)
{
    std::random_device dev_rnd;
    uint64_t initstate = dev_rnd();
    uint64_t initseq = dev_rnd();
    pcg32_srandom_r(rng, initstate, initseq);
}

float rnd(pcg32_random_t* rng) //return [0,1)
{
    return std::ldexp(pcg32_random_r(rng), -32);
}

#endif