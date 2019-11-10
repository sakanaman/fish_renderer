#include "concurrent.hpp"

ParallelRender::ParallelRender(const std::function<void(int*, int*)>& _render):render(_render), count(0){}


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
        [&]()
        {
            int upper_left[2] = {-1.0, -1.0};
            int bottom_right[2] = {-1.0, -1.0};

            while(1)
            {
                {
                    std::lock_guard<std::mutex> lg(mut);
                    if(count < split_num * split_num)
                    {
                        ++count;
                        int i = ((count % split_num) + split_num - 1) % split_num;
                        int j = count / split_num;
                        upper_left[0] = w[i];
                        upper_left[1] = h[j];
                        bottom_right[0] = w[i + 1];
                        bottom_right[1] = h[j + 1];

                    }
                    else
                    {
                        return;
                    }
                }

                render(upper_left, bottom_right);
            }
        };

    int thread_num = std::thread::hardware_concurrency();

    std::vector<std::thread> threads(thread_num);
    for(int i = 0; i < thread_num; ++i)
    {
        threads[i] = std::thread(renderthread);
    }

    for(int i = 0; i < thread_num; i++)
    {
        threads[i].join();
    }

    std::cout << count << std::endl;
}