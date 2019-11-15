#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <Vec.hpp>
#include <random.hpp>

template<class Real>
class MaterialData
{
public:
    MaterialData(){};
    MaterialData(const Vec3<Real>& speculer, const Vec3<Real>& diffuse, 
                  const  Vec3<Real>& glossy, const Vec3<Real>& transmit, const Vec3<Real> emission,
                  Real roughness, Real ni):
                  speculer(speculer), diffuse(diffuse), transmit(transmit),
                  emission(emission),roughness(roughness), ni(ni){}
    Vec3<Real> GetSpeculer()
    {  
        return speculer;
    };
    Vec3<Real> GetDiffuse()
    {
        return diffuse;
    };
    Vec3<Real> GetTrasimit()
    {
        return transmit;
    };
    Real GetRoughness()
    {
        return roughness;
    };
    Real GetNi()
    {
        return ni;
    };
private:
    Vec3<Real> speculer;
    Vec3<Real> diffuse;
    Vec3<Real> transmit;
    Vec3<Real> emission;
    Real roughness;
    Real ni;
};

// template<class Real>
// Vec3<Real> fresnel(const Vec3<Real>& speculer, const Vec3<Real>& wo, const Vec3<Real>& n)
// {
//     Real costerm = std::abs(wo.dot(n));
//     return speculer + (Vec3<Real>(1.0,1.0,1.0) - speculer) * std::pow(1 - costerm, 5.0f);
// }

// template<class Real>
// Real SmithG1(const Vec3<Real>& v, const Vec3<Real>& m, Real alpha){
//     Real costhetaI = v.dot(m);
//     Real sinthetaI = std::sqrt(1 - costhetaI * costhetaI);
//     Real tanthetaI = sinthetaI/costhetaI;
//     return 2.0/(1.0 + std::sqrt(1.0 + alpha*alpha*tanthetaI*tanthetaI));
// }

// template<class Real>
// Real SmithG2(const Vec3<Real>& wi, const Vec3<Real>& wo, const Vec3<Real>& n, Real alpha)
// {
//     Real Dotgo = wo.dot(n);
//     Real Dotgi = wi.dot(n);
//     Real G2_denom1 = Dotgo*std::sqrt(alpha * alpha + (1 - alpha*alpha)*Dotgi*Dotgi);
//     Real G2_denom2 = Dotgi*std::sqrt(alpha * alpha + (1 - alpha*alpha)*Dotgo*Dotgo);
    
//     return 2.0 * Dotgi * Dotgo / (G2_denom1 + G2_denom2);
// }

// template<class Real>
// Vec3<Real> SampleGGX(const Vec3<Real>& wo, const MaterialData<Real>& mat , const Vec3<Real>& n, 
//                     Vec3<Real>& weight, pcg32_random_t* rng_ptr)
// {
//     Real roughness = mat.GetRoughness();
//     Real alpha = roughness * roughness;

//     Vec3<Real> e0, e2;
//     ONB(n, e0, e2);
//     Vec3<Real> local_wi = worldtolocal(wo, e0, n, e2);

//     Vec3<Real> wh, wi;

//     do {
//         Vec3<Real> V = Vec3<Real>(alpha*local_wi[0], local_wi[1], alpha * local_wi[2]).normalize();

//         Vec3<Real> T1 =  (V[1] < 0.9999) ? (V.cross(Vec3<Real>(0,1,0))).normalize() : Vec3<Real>(1,0,0);
// 		Vec3<Real> T2 = T1.cross(V);

//         Real U1 = static_cast<Real>(rnd(rng_ptr)), U2 = static_cast<Real>(rnd(rng_ptr));

//         Real a = 1.0 / (1.0 + V[1]);
// 		Real r = sqrt(U1);
// 		Real phi = (U2<a) ? U2/a * M_PI : M_PI + (U2-a)/(1.0-a) * M_PI;
// 		Real P1 = r*cos(phi);
// 		Real P2 = r*sin(phi)*((U2<a) ? 1.0 : V[1]);

//         Vec3<Real> N = P1*T1 + P2*T2 + sqrt(std::max(0.0, 1.0 - P1*P1 - P2*P2))*V;
//         N = (Vec3<Real>(alpha*N[0], std::max(0.0, N[1]), alpha*N[2])).normalize();

//         wh = (N[0] * e0 + N[1] * n + N[2] * e2).normalize();
//     }while( wh.dot(n) < 0);

//     wi = 2.0 * std::abs(wo.dot(wh))*wh - wo;

//     Vec3<Real> speculer = mat.GetSpeculer();
//     Vec3<Real> fr = fresnel(speculer, wi, wh);
//     weight = fr * SmithG2(wi, wo, n, alpha) / SmithG1(wo, n, alpha);

//     return wi;
// }

// template<class Real>
// Vec3<Real> SampleDiffuse(const Vec3<Real>& wo, const MaterialData<Real>& mat , const Vec3<Real>& n, 
//                         Vec3<Real>& weight, pcg32_random_t* rng_ptr)
// {
//     weight = mat.GetDiffuse();

//     //sample
//     Vec3<Real> e0, e2;
//     ONB(n, e0, e2);
//     Real u1 = rnd(rng_ptr), u2 = rnd(rng_ptr);

//     Real y = u1;
//     Real x = std::sqrt(1 - y*y) * std::cos(2*M_PI*u2);
//     Real z = std::sqrt(1 - y*y) * std::sin(2*M_PI*u2);
//     return x * e0 + y * n + z * e2;
// }

// //layer material
// template<class Real>
// Vec3<Real> SampleGrosDiffuse(const Vec3<Real>& wo, const MaterialData<Real>& mat , const Vec3<Real>& n, 
//                              Vec3<Real>& weight, pcg32_random_t* rng_ptr)  //Grossy + Diffuse
// {
//     Vec3<Real> speculer = mat.GetSpeculer();
//     Vec3<Real> fr = fresnel(speculer, wo, n);
//     Real p1 = static_cast<Real>(rnd(rng_ptr));
    
//     Vec3<Real> wi;
//     if(p1 < (fr[0] + fr[1] + fr[2])/3.0) //Glossy(GGX)
//     {
//         wi = SampleDiffuse(wo, mat, n, weight, rng_ptr);   
//     }
//     else // Diffuse
//     {   
//         wi = SampleGGX(wo, mat, n, weight, rng_ptr);
//     }
//     return wi;
// } //this function return "wi"

template class MaterialData<float>;
#endif