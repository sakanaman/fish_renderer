#ifndef RAYTRACE_HPP
#define RAYTRACE_HPP

#include "material.hpp"

class Ray
{
public:
    Ray(){}
    Ray(const Vec3<float>& origin, const Vec3<float>& dir)
    :origin(origin),dir(dir){}
    Ray(const float* _origin, const float* _dir)
    {
        origin = Vec3<float>(_origin);
        dir = Vec3<float>(_dir);
    }

    Vec3<float> origin;
    Vec3<float> dir;
};

Vec3<float> Trace_Test(const float* firstRay_dir, const float* firstRay_origin, const Shader<float>& shader, const IBL& ibl,
                       const BinaryBVH<float, TriangleData<float>>& bvh, RandomManager& rnd_manager)
{
    Vec3<float> throughput(1.0f, 1.0f, 1.0f);
    Vec3<float> I(0.0f, 0.0f, 0.0f);
    Ray comingRay(firstRay_origin, firstRay_dir);
    float Pr = 1.0;
    for(int depth = 0;;++depth)
    {
        Hit<float> hit;
        float comingRay_origin[] = {comingRay.origin[0], comingRay.origin[1], comingRay.origin[2]};
        float comingRay_dir[] = {comingRay.dir[0], comingRay.dir[1], comingRay.dir[2]};
        bool is_hit = bvh.Traverse(hit, comingRay_origin, comingRay_dir, 0);
        
        if(is_hit)
        {
            Vec3<float> wo = -1.0f * comingRay.dir;
            Vec3<float> n(hit.GetNg()[0], hit.GetNg()[1], hit.GetNg()[2]);
            Vec3<float> hitPos(hit.GetPos()[0], hit.GetPos()[1], hit.GetPos()[2]);
            Vec3<float> orienting_normal = n.dot(comingRay.dir) < 0 ? n : -1.0f * n;

            Vec3<float> e0, e2;
            ONB(orienting_normal, e0, e2);
            Vec3<float> local_wo{e0.dot(wo),
                                 orienting_normal.dot(wo),
                                 e2.dot(wo)};
            

            #ifdef NEE
            float phi, theta;
            float Le[3] = {0.0f, 0.0f, 0.0f};
            float r1 = rnd_manager.GetRND(), r2 = rnd_manager.GetRND();
            float nee_pdf = ibl.sample(r1, r2, &phi, &theta, Le);

            Vec3<float> shadowdir = {std::sin(theta) * std::cos(phi),
                                     std::cos(theta),
                                     std::sin(theta) * std::sin(phi)};
            Vec3<float> local_shadowdir = {e0.dot(shadowdir),
                                           orienting_normal.dot(shadowdir),
                                           e2.dot(shadowdir)};
            
            if(shadowdir.dot(orienting_normal) > 0.0f)
            {
                float shadowdir_ary[3] = {shadowdir[0], shadowdir[1], shadowdir[2]};
                Vec3<float> shadoworigin  = hitPos /*+ 0.001f * orienting_normal*/;
                float shadoworigin_ary[3] = {shadoworigin[0], shadoworigin[1], shadoworigin[2]};
                Hit<float> shadow_hit;
                if(bvh.Traverse(shadow_hit, shadoworigin_ary, shadowdir_ary, 0))
                {
                    //nothing
                }
                else
                {
                    float costerm = local_shadowdir[1];
                    I = I + throughput * eval_bsdf(local_shadowdir, local_wo, hit, shader) 
                                       * (1.0f/nee_pdf) * costerm * Vec3<float>(Le);
                }
            }
            #endif

            Vec3<float> local_wi = sample(local_wo, hit, rnd_manager, shader);
            Vec3<float> wi = local_wi[0]*e0 + 
                             local_wi[1]*orienting_normal + 
                             local_wi[2]*e2;
            comingRay = Ray(hitPos /*+ 0.001f * orienting_normal*/, wi);

            throughput = throughput * weight(local_wi, local_wo, hit, shader);

            Pr *= 0.96;
            if(rnd_manager.GetRND() < Pr)
            {
                throughput  = throughput * (1.0f/Pr);
            }
            else
            {
                break;
            }
        }
        else
        {
            #ifdef NEE
            if(depth == 0)
            {
                float phi = std::atan2(comingRay.dir[2], comingRay.dir[0]);
                float theta = std::acos(comingRay.dir[1]);
                if (phi < 0)
                    phi += 2 * M_PI;
                if (phi > 2 * M_PI)
                    phi -= 2 * M_PI;

                float Le[3];
                ibl.GetLe(Le, theta, phi);
                I = I + throughput * Vec3<float>(Le);
            }
            #else
                float phi = std::atan2(comingRay.dir[2], comingRay.dir[0]);
                float theta = std::acos(comingRay.dir[1]);
                if (phi < 0)
                    phi += 2 * M_PI;
                if (phi > 2 * M_PI)
                    phi -= 2 * M_PI;

                float Le[3];
                ibl.GetLe(Le, theta, phi);
                I = I + throughput * Vec3<float>(Le);
            #endif
            break;
        }   
    }
    return I;
}

Vec3<float> Trace(const float* firstRay_dir, const float* firstRay_origin, const SceneData<float>& scenedata, const IBL& ibl,
                  const BinaryBVH<float, TriangleData<float>>& bvh, RandomManager& rnd_manager)
{
    Vec3<float> throughput(1.0f, 1.0f, 1.0f);
    Vec3<float> I(0.0f, 0.0f, 0.0f);
    Ray comingRay(firstRay_origin, firstRay_dir);
    float Pr = 1.0;
    for(int depth = 0;;++depth)
    {
        Hit<float> hit;
        float comingRay_origin[] = {comingRay.origin[0], comingRay.origin[1], comingRay.origin[2]};
        float comingRay_dir[] = {comingRay.dir[0], comingRay.dir[1], comingRay.dir[2]};
        bool is_hit = bvh.Traverse(hit, comingRay_origin, comingRay_dir, 0);
        
        if(is_hit)
        {
            Vec3<float> n(hit.GetNg()[0], hit.GetNg()[1], hit.GetNg()[2]);
            Vec3<float> hitPos(hit.GetPos()[0], hit.GetPos()[1], hit.GetPos()[2]);
            Vec3<float> orienting_normal = n.dot(comingRay.dir) < 0 ? n : -1.0f * n;

            Vec3<float> e0, e2;
            ONB(orienting_normal, e0, e2);

            //Light Sampling
            #ifdef NEE
            float phi, theta;
            float Le[3] = {0.0f, 0.0f, 0.0f};
            float r1 = rnd_manager.GetRND(), r2 = rnd_manager.GetRND();
            float nee_pdf = ibl.sample(r1, r2, &phi, &theta, Le);
            
            //std::cout << Vec3<float>(Le) << std::endl;
            Vec3<float> shadowdir = std::sin(theta) * std::cos(phi) * Vec3<float>(1.0f, 0.0f, 0.0f)
                                  + std::cos(theta) * Vec3<float>(0.0f, 1.0f, 0.0f)
                                  + std::sin(theta) * std::sin(phi) * Vec3<float>(0.0f, 0.0f, 1.0f);
            if(shadowdir.dot(orienting_normal) > 0.0f)
            {
                float shadowdir_ary[3] = {shadowdir[0], shadowdir[1], shadowdir[2]};
                Vec3<float> shadoworigin  = hitPos /*+ 0.001f * orienting_normal*/;
                float shadoworigin_ary[3] = {shadoworigin[0], shadoworigin[1], shadoworigin[2]};
                Hit<float> shadow_hit;
                if(bvh.Traverse(shadow_hit, shadoworigin_ary, shadowdir_ary, 0))
                {
                    //nothing
                }
                else
                {
                    float costerm = orienting_normal.dot(shadowdir);
                    I[0] = I[0] + throughput[0] * scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[hit.GetID()] + 0] * Le[0]
                                                * costerm
                                                / (nee_pdf * M_PI);
                    I[1] = I[1] + throughput[1] * scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[hit.GetID()] + 1] * Le[1]
                                                * costerm
                                                / (nee_pdf * M_PI);
                    I[2] = I[2] + throughput[2] * scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[hit.GetID()] + 2] * Le[2]
                                                * costerm
                                                / (nee_pdf * M_PI);
                }
            }
            #endif

            //BRDF sampling
            float u1 = rnd_manager.GetRND();
            float u2 = rnd_manager.GetRND();
            Vec3<float> wi = randomCosineHemisphere(u1, u2, orienting_normal);

            comingRay = Ray(hitPos /*+ 0.001f * orienting_normal*/, wi);

            // store throughput
            throughput[0] *= scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[hit.GetID()] + 0];
            throughput[1] *= scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[hit.GetID()] + 1];
            throughput[2] *= scenedata.mat_infos.diffuses[3 * scenedata.mat_infos.mat_indices[hit.GetID()] + 2];

            // float Pr = std::max({throughput[0], throughput[1], throughput[2]});
            //russian roulette
            Pr *= 0.96;
            if(rnd_manager.GetRND() < Pr)
            {
                throughput  = throughput * (1.0f/Pr);
            }
            else
            {
                break;
            }
            
        }
        else
        {
            #ifdef NEE
            if(depth == 0)
            {
                float phi = std::atan2(comingRay.dir[2], comingRay.dir[0]);
                float theta = std::acos(comingRay.dir[1]);
                if (phi < 0)
                    phi += 2 * M_PI;
                if (phi > 2 * M_PI)
                    phi -= 2 * M_PI;

                float Le[3];
                ibl.GetLe(Le, theta, phi);
                I = I + throughput * Vec3<float>(Le);
            }
            #else
                float phi = std::atan2(comingRay.dir[2], comingRay.dir[0]);
                float theta = std::acos(comingRay.dir[1]);
                if (phi < 0)
                    phi += 2 * M_PI;
                if (phi > 2 * M_PI)
                    phi -= 2 * M_PI;

                float Le[3];
                ibl.GetLe(Le, theta, phi);
                I = I + throughput * Vec3<float>(Le);
            #endif
            break;
        }   
    }
    return I;
}


Vec3<float> Trace_debug(const float* firstRay_dir, const float* firstRay_origin, const SceneData<float>& scenedata, const IBL& ibl,
                        const BinaryBVH<float, TriangleData<float>>& bvh, RandomManager& rnd_manager)
{
    Vec3<float> light_dir = Vec3<float>(-1.0f, 2.0f, 0.0f).normalize();

    Vec3<float> throughput(1.0f, 1.0f, 1.0f);
    Vec3<float> I(0.0f, 0.0f, 0.0f);
    Ray comingRay(firstRay_origin, firstRay_dir);

    Hit<float> hit;
    float comingRay_origin[] = {comingRay.origin[0], comingRay.origin[1], comingRay.origin[2]};
    float comingRay_dir[] = {comingRay.dir[0], comingRay.dir[1], comingRay.dir[2]};
    bool is_hit = bvh.Traverse(hit, comingRay_origin, comingRay_dir, 0);

    if(is_hit)
    {
        Vec3<float> n(hit.GetNg()[0], hit.GetNg()[1], hit.GetNg()[2]);
        Vec3<float> hitPos(hit.GetPos()[0], hit.GetPos()[1], hit.GetPos()[2]);
        Vec3<float> orienting_normal = n.dot(comingRay.dir) < 0 ? n : -1.0f * n;
        Vec3<float> norm_color{(n[0] + 1)*0.5f,
                               (n[1] + 1)*0.5f,
                               (n[2] + 1)*0.5f};


        Vec3<float> shadowdir = light_dir;
        float shadowdir_ary[3] = {shadowdir[0], shadowdir[1], shadowdir[2]};
        Vec3<float> shadoworigin  =  hitPos;
        float shadoworigin_ary[3] = {shadoworigin[0], shadoworigin[1], shadoworigin[2]};
        Hit<float> shadow_hit;
        if(bvh.Traverse(shadow_hit, shadoworigin_ary, shadowdir_ary, 0))
        {
            if(shadow_hit.GetID() == hit.GetID()) //oh..self intersection...
            {
                return {0.0f, 0.0f, 0.0f};
            }
            return {0.3f * norm_color[0], 
                    0.3f * norm_color[1],
                    0.3f * norm_color[2]};
        }
        return {norm_color[0], 
                norm_color[1],
                norm_color[2]};
    }
    else
    {
        return {0.0f, 0.0f, 0.0f};
    }
}

#endif