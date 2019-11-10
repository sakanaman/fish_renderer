#ifndef CONCURRENT_HPP
#define CONCURRENT_HPP

#include <vector>
#include <functional>
#include <mutex>
#include <atomic>
#include <thread>
#include <iostream>

//task allocate with tile splitt
class ParallelRender
{
public:
    ParallelRender(const std::function<void(int*, int*)>& _render);
    void Execute(const int width, const int height, const int split_num);
private:
    std::atomic<int> count;
    std::mutex mut;
    std::function<void(int*, int*)> render;
    // --> render function argument = (upper_left[2], bottom_right[2])
};

#endif