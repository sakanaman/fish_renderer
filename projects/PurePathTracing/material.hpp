#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <variant>
#include "scene.hpp"
#include "bsdf.hpp"
#include "sampling.hpp"

enum class ShaderType
{
    R, // reflection only
    TR, // Transmit or Reflection
    TRS, // Reflection or Transmit & SSS (Mainly for Uber Shader)
    VOL, // volume only
};

//diffuse only shader(for test)
template<class Real>
class DiffuseShader
{
public:
    DiffuseShader(const SceneData<Real>& scenedata)
                  :scenedata(scenedata){};
    //for surface
    Vec3<Real> sample(const Vec3<Real>& wo, const Hit<Real>& hit, RandomManager& rnd_manager) const
    {
        Real u = rnd_manager.GetRND();
        Real v = rnd_manager.GetRND();
        return test_randomCosineHemisphere(u, v);
    };
    Real eval_pdf(const Vec3<Real>& wi, const Vec3<Real>& wo,const Hit<Real>& hit) const
    {
        return eval_diffuse_pdf(wi, wo);
    };
    Vec3<Real> eval_bsdf(const Vec3<Real>& wi, const Vec3<Real>& wo, const Hit<Real>& hit) const
    {
        int geomID = hit.GetID();
        Vec3<Real> albedo{scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 0],
                          scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 1],
                          scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 2]};
        return eval_diffuse_bsdf(wi, wo, albedo);
    };
    Vec3<Real> weight(const Vec3<Real>& wi, const Vec3<Real>& wo, const Hit<Real>& hit) const
    {
        int geomID = hit.GetID();
        Vec3<Real> albedo{scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 0],
                          scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 1],
                          scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[geomID] + 2]};
        return albedo;
    };

    //for subsurface or vpt
    Vec3<Real> sample_phase(const Vec3<Real>& wo, const Hit<Real>& hit){};
    void CalcDensity(const Vec3<Real>& position, Real* sigma_t, Real* sigma_s){}; //this shader is not for vpt
private:
    SceneData<Real> scenedata;
    const ShaderType shader_type = ShaderType::R;
};


//unimplemented
// template<class Real>
// class SimpleUberShader
// {
// public:
//     //for surface
//     Vec3<Real> sample(const Vec3<Real>& wo, const Hit<Real>& hit, RandomManager& rnd_manager){};
//     Real eval_pdf(const Vec3<Real>& wi, const Vec3<Real>& wo,const Hit<Real>& hit){};
//     Vec3<Real> eval_bsdf(const Vec3<Real>& wi, const Vec3<Real>& wo, const Hit<Real>& hit){};
//     Vec3<Real> weight(const Vec3<Real>& wi, const Vec3<Real>& wo, const Hit<Real>& hit){};

//     //for subsurface or vpt
//     Vec3<Real> sample_phase(const Vec3<Real>& wo, const Hit<Real>& hit){};
//     void CalcDensity(const Vec3<Real>& position, Real* sigma_t, Real* sigma_s){}; //this shader is not for vpt
// private:
//     SceneData<Real> scenedata;
//     const ShaderType shader_type = ShaderType::TR;
// };

// //unimplemented
// template<class Real>
// class Disney2012Shader
// {
// public:
//     //for surface
//     Vec3<Real> sample(const Vec3<Real>& wo, const Hit<Real>& hit,RandomManager& rnd_manager){};
//     Real eval_pdf(const Vec3<Real>& wi, const Vec3<Real>& wo,const Hit<Real>& hit){};
//     Vec3<Real> eval_bsdf(const Vec3<Real>& wi, const Vec3<Real>& wo, const Hit<Real>& hit){};
//     Vec3<Real> weight(const Vec3<Real>& wi, const Vec3<Real>& wo, const Hit<Real>& hit){};

//     //for subsurface or vpt
//     Vec3<Real> sample_phase(const Vec3<Real>& wo, const Hit<Real>& hit){};
//     void CalcDensity(const Vec3<Real>& position, Real* sigma_t, Real* sigma_s){}; //this shader is not for vpt
// private:
//     SceneData<Real> scenedata;
//     const ShaderType shader_type = ShaderType::R;
// };


template<class T>
using Shader = std::variant<DiffuseShader<T>>;


// wrapper functions 
template<class Real>
Vec3<Real> sample(const Vec3<Real>& wo, const Hit<Real>& hit, RandomManager& rnd_manager
                 ,const Shader<Real>& shader)
{
    return std::visit([&](auto& x){return x.sample(wo, hit, rnd_manager);}, shader);
}

template<class Real>
Real eval_pdf(const Vec3<Real>& wi, const Vec3<Real>& wo,const Hit<Real>& hit
             ,const Shader<Real>& shader)
{
    return std::visit([&](auto& x){return x.eval_pdf(wi, wo, hit);},shader);
}

template<class Real>
Vec3<Real> eval_bsdf(const Vec3<Real>& wi, const Vec3<Real>& wo, const Hit<Real>& hit
                    ,const Shader<Real>& shader)
{
    return std::visit([&](auto& x){return x.eval_bsdf(wi, wo, hit);}, shader);
}

template<class Real>
Vec3<Real> weight(const Vec3<Real>& wi, const Vec3<Real>& wo, const Hit<Real>& hit,
                  const Shader<Real>& shader)
{
    return std::visit([&](auto& x){return x.weight(wi, wo, hit);}, shader);
}
#endif