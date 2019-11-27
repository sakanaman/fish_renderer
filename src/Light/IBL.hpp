#ifndef IBL_HPP
#define IBL_HPP
#include <stb_image.h>
#include <vector>
#include <string>


class IBL
{
public:
    IBL(std::string _filename);
    void xy2uv(const int x, const int y, float* u, float* v) const;
    void uv2xy(int* x, int* y, const float u, const float v) const;
    void uv2angle(const float u, const float v, float* theta, float* phi) const;
    void angle2uv(float* u, float* v, const float theta, const float phi) const;
    void GetLe(float* Le, const float theta, const float phi) const;
    float sample(const float r1, const float r2, float* phi, float* theta, float* Le) const;
private:
    std::string name;
    float* data;
    int channel;
    int w;
    int h;

    //for importance sampling
    std::vector<std::vector<float>> PVDists;
    std::vector<float> PUDists;
};

#endif