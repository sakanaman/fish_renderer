#ifndef CONCURRENT_HPP
#define CONCURRENT_HPP

#include <vector>
#include <functional>
#include <mutex>
#include <atomic>
#include <thread>
#include <iostream>
#include "random.hpp"

//task allocate with tile splitt
class ParallelRender
{
public:
    ParallelRender(const std::function<void(const int*,const int*, pcg32_random_t*)>& _render);
    void Execute(const int width, const int height, const int split_num);
private:
    std::atomic<int> count;
    std::mutex mut;
    std::function<void(const int*, const int*, pcg32_random_t*)> render;
    // --> render function argument = (upper_left[2], bottom_right[2], rng_ptr)
};

//
ParallelRender::ParallelRender(const std::function<void(const int*, const int*, pcg32_random_t*)>& _render):count(0),render(_render){}


void ParallelRender::Execute(const int width, const int height, const int split_num)
{
    // get split point: [ w[i], w[i + 1] ) × [ h[j], h[j + 1] ) = a task  (i,j = 0,1,2,...split_num-1)
    int w[split_num + 1], h[split_num + 1];
    int w_chunk = width/split_num, h_chunk = height/split_num;
    for(int i = 0; i < split_num; ++i)
    {
        w[i] = i * w_chunk;
        h[i] = i * h_chunk;
    }
    h[split_num] = height;
    w[split_num] = width;

    //     i    
    //  ________
    //  |
    // j|
    //  |
    auto renderthread = 
        [&](pcg32_random_t* rng)
        {
            int upper_left[2] = {-1, -1};
            int bottom_right[2] = {-1, -1};

            while(1)
            {
                bool have_task = false;
                { // task allocating must be locked !!
                    std::lock_guard<std::mutex> lg(mut);
                    if(count < split_num * split_num)
                    {
                        int i = ((count % split_num) + split_num - 1) % split_num;
                        int j = count / split_num;
                        upper_left[0] = w[i];
                        upper_left[1] = h[j];
                        bottom_right[0] = w[i + 1];
                        bottom_right[1] = h[j + 1];
                        ++count;
                        fprintf(stderr, "\r[%3d]", int(100 * count / (split_num * split_num))); 
                        have_task = true;
                    }
                }
                if(!have_task) return;

                render(upper_left, bottom_right, rng);
            }
        };

    int thread_num = std::thread::hardware_concurrency();
    std::vector<pcg32_random_t> rngs(thread_num);
    for(int i = 0; i < thread_num; ++i)
    {
        initRNG(&(rngs[i]));
    }

    std::vector<std::thread> threads(thread_num);
    for(int i = 0; i < thread_num; ++i)
    {
        threads[i] = std::thread(renderthread, &(rngs[i]));
    }

    for(int i = 0; i < thread_num; i++)
    {
        threads[i].join();
    }

    std::cout << count << std::endl;
}
//

#endif