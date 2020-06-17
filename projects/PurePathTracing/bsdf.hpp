#ifndef BSDF_HPP
#define BSDF_HPP

#include "api.hpp"
#include <cmath>

//Diffuse
template<class Real>
Vec3<Real> eval_diffuse_bsdf(const Vec3<Real>& wi, const Vec3<Real>& wo, 
                             Vec3<Real> albedo)
{
    return albedo * Real(1.0/M_PI);
}

template<class Real>
Real eval_diffuse_pdf(const Vec3<Real>& wi, const Vec3<Real>& wo)
{
    Real costheta = wi.y; //because y-up local coord
    return costheta * (1.0/M_PI);
}

template<class Real>
Vec3<Real> eval_diffuse_weight(const Vec3<Real>& wi, const Vec3<Real>& wo, 
                             Vec3<Real> albedo)
{
    return albedo; 
}                            


#endif