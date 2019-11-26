#ifndef CONCURRENT_HPP
#define CONCURRENT_HPP

#include <vector>
#include <functional>
#include <mutex>
#include <atomic>
#include <thread>
#include <iostream>
#include "random.hpp"

class RandomManager
{
public:
    RandomManager();
    float GetRND();
private:
    pcg32_random_t rnd_class;
};

//task allocate with tile splitt
class ParallelRender
{
public:
    ParallelRender(const std::function<void(const int*,const int*, RandomManager&)>& _render);
    void Execute(const int width, const int height, const int split_num);
private:
    std::atomic<int> count;
    std::mutex mut;
    std::function<void(const int*, const int*, RandomManager&)> render;
    // --> render function argument = (upper_left[2], bottom_right[2], rng_ptr)
};

//
//

#endif