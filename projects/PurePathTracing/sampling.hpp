#ifndef SAMPLING_HPP
#define SAMPLING_HPP

#include "api.hpp"

// Diffuse
Vec3<float> randomCosineHemisphere(const float u, const float v, const Vec3<float> &n)
{
    float theta = 0.5 * std::acos(1 - 2 * u);
    float phi = 2 * M_PI * v;

    float x = std::cos(phi) * std::sin(theta);
    float y = std::cos(theta);
    float z = std::sin(phi) * std::sin(theta);
    Vec3<float> xv, zv;
    ONB(n, xv, zv);
    return x * xv + y * n + z * zv;
}

Vec3<float> test_randomCosineHemisphere(const float u, const float v)
{
    float theta = 0.5 * std::acos(1 - 2 * u);
    float phi = 2 * M_PI * v;

    float x = std::cos(phi) * std::sin(theta);
    float y = std::cos(theta);
    float z = std::sin(phi) * std::sin(theta);
    return {x, y, z};
}

#endif