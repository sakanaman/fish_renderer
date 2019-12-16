#include "bssrdf.hpp"
#include <cmath>


//////Material Infomation Zone//////
float FresnelMoment1(float eta) 
{
    float eta2 = eta * eta, eta3 = eta2 * eta, eta4 = eta3 * eta,
          eta5 = eta4 * eta;
    if (eta < 1)
        return 0.45966f - 1.73965f * eta + 3.37668f * eta2 - 3.904945 * eta3 +
               2.49277f * eta4 - 0.68441f * eta5;
    else
        return -4.61686f + 11.1136f * eta - 10.4646f * eta2 + 5.11455f * eta3 -
               1.27198f * eta4 + 0.12746f * eta5;
}

float Sw(float eta, float costerm)
{
    float c = 1 - 2 * FresnelMoment1(1 / eta);
    float R = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
    float fres_diel = R + (1.0f - R) * std::pow(1.0f - costerm, 5.0);
    return (1 - fres_diel) / (c * M_PI);
}

Vec3<float> CalcScallFactor(const Vec3<float>& albedo)
{
    Vec3<float> s;
    s[0] = 1.85f - albedo[0] + 7.0f * std::pow(std::abs(albedo[0] - 0.8f),3.0f);
    s[1] = 1.85f - albedo[1] + 7.0f * std::pow(std::abs(albedo[1] - 0.8f),3.0f);
    s[2] = 1.85f - albedo[2] + 7.0f * std::pow(std::abs(albedo[2] - 0.8f),3.0f);
    return s;
}

Vec3<float> CalcCurveParam(const Vec3<float>& s)
{
    return Vec3<float>(1.0f, 1.0f, 1.0f) / s; // this is ok! since mfp=1
}

Vec3<float> Sr(float r, const Vec3<float>& albedo, const Vec3<float>& s)
{
    if(r < 1e-6f) r = 1e-6f;
     
    Vec3<float> R;
    R[0] = albedo[0] * s[0] * (std::exp(-s[0]*r) + std::exp(-s[0]*r/3.0f))
                            / (8.0f * M_PI * r);
    R[1] = albedo[1] * s[1] * (std::exp(-s[1]*r) + std::exp(-s[1]*r/3.0f))
                            / (8.0f * M_PI * r);
    R[2] = albedo[2] * s[2] * (std::exp(-s[2]*r) + std::exp(-s[2]*r/3.0f))
                            / (8.0f * M_PI * r);

    return R;
}

float Sample_Sr(int ch, float u, const Vec3<float>& albedo, const Vec3<float>& d)
{
    const float OneMinusEpsilon = 0.99999994;
    if(u < 0.25f)
    {
        u = std::min<float>(u * 4, OneMinusEpsilon);
        return d[ch] * std::log(1/(1-u));
    }
    else
    {
        u = std::min<float>((u - .25f) / .75f, OneMinusEpsilon);
        return 3 * d[ch] * std::log(1/(1-u));
    }
}

float Pdf_Sr(int ch, float r, const Vec3<float>& d)
{
    if(r < 1e-6f) r = 1e-6f;

    return (.25f * std::exp(-r / d[ch]) / (2 * M_PI * d[ch] * r) +
            .75f * std::exp(-r / (3 * d[ch])) / (6 * M_PI * d[ch] * r));
}


////// Implementation Zone ///////
// Axis MIS -> Sample_Sr
Vec3<float> Sample_Sp(const Vec3<float>& p0, const Vec3<float>& wo, // p0 data
                      const Vec3<float>& Albedo,
                      const BinaryBVH<float, TriangleData<float>>& bvh, // scene data
                      const float u1, const float u2, const float u3, const float u4, const float u5,// random number's !
                      const Vec3<float>& T, const Vec3<float>& N, const Vec3<float> B, // orthonormal-basis
                      Hit<float>* Probehit, float* pdf)
{
    // Axis Sampling
    Vec3<float> vx, vy, vz;
    if(u1 < 0.5f)
    {
        vx = T;
        vy = B;
        vz = N;
    }
    else if(u1 < 0.75f)
    {
        vx = B;
        vy = N;
        vz = T;
    }
    else
    {
        vx = N;
        vy = T;
        vz = B;
    }

    // Select color channel from rgb
    int ch = u2 * 3;

    // sample r
    Vec3<float> s = CalcScallFactor(Albedo);
    Vec3<float> d = CalcCurveParam(s); 
    float r = Sample_Sr(ch, u3, Albedo, d);
    if(r < 0)
    {
        return Vec3<float>(0.0f, 0.0f, 0.0f);
    }
    float phi = 2 * M_PI * u4;

    float rMax = Sample_Sr(ch, 0.999f, Albedo, d);
    if(r >= rMax) return Vec3<float>(0.0f, 0.0f, 0.0f);
    float l = 2 * std::sqrt(rMax * rMax - r * r);

    // track ProbeRay
    std::vector<Hit<float>> ProbeHits;
    Vec3<float> start = p0 + vx * r * cos(phi) + vy * r * sin(phi) - l * 0.5f * vz;
    float dir[] = {vz[0], vz[1], vz[2]};
    Vec3<float> modify_start;
    float accum_dist = 0.0f;
    while(true)
    {   
        Hit<float> Phit;
        modify_start = start + 0.01f * vz;
        float o[] = {modify_start[0], modify_start[1], modify_start[2]};
        bool is_probehit = bvh.Traverse(Phit, o, dir, 0);
        if(is_probehit && (accum_dist + Phit.GetT() < l))
        {
            ProbeHits.push_back(Phit);
            start = Vec3<float>(Phit.GetPos());
            accum_dist += Phit.GetT();
        }
        else
        {
            break;
        }
    }
    if(ProbeHits.empty()) return Vec3<float>(0.0f, 0.0f, 0.0f);

    // select incident point
    int winner_p1 = u5 * ProbeHits.size();
    *Probehit = ProbeHits[winner_p1];

    //Pdf
    Vec3<float> n_i(Probehit->GetNg());
    float pdf_axis[3] = {0.5f, 0.25f, 0.25f};
    float axis_cos[3] = {n_i.dot(N), n_i.dot(T), n_i.dot(B)};
    float pdf_rgb[3] = {1.0f/3.0f, 1.0f/3.0f, 1.0f/3.0f};
    for(int color_ch = 0; color_ch < 3; color_ch++)
    {
        for(int axis = 0; axis < 3; axis++)
        {
            *pdf += pdf_rgb[color_ch] * pdf_axis[axis] * Pdf_Sr(color_ch, r, d) * std::abs(axis_cos[axis]);
                    //rgb               //axis           //disk                   //cos
        }
    }
    *pdf /= ProbeHits.size();

    //Sr(diffusion term)
    Vec3<float> sr = Sr(r, Albedo, s);
    return sr;
}