#ifndef BSSRDF_HPP
#define BSSRDF_HPP
#include <Vec.hpp>
#include <Accelerator.hpp>
float FresnelMoment1(float eta);
float Sw(float eta, float costerm);
Vec3<float> CalcScallFactor(const Vec3<float>& albedo);
Vec3<float> CalcCurveParam(const Vec3<float>& s);
Vec3<float> Sr(float r, const Vec3<float>& albedo, const Vec3<float>& s);
float Sample_Sr(int ch, float u, const Vec3<float>& albedo, const Vec3<float>& d);
float Pdf_Sr(int ch, float r, const Vec3<float>& d);
Vec3<float> Sample_Sp(const Vec3<float>& p0, const Vec3<float>& wo, // p0 data
                      const Vec3<float>& Albedo,
                      const BinaryBVH<float, TriangleData<float>>& bvh, // scene data
                      const float u1, const float u2, const float u3, const float u4, const float u5,// random number's !
                      const Vec3<float>& T, const Vec3<float>& N, const Vec3<float> B, // orthonormal-basis
                      Hit<float>* Probehit, float* pdf);

#endif