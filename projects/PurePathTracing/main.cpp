#include <camera.hpp>
#include <Accelerator.hpp>
#include <tiny_obj_loader.h>
#include <stb_image_write.h>
#include <IBL.hpp>
#include <concurrent.hpp>
#include <Vec.hpp>
#define NEE

template<class Real>
class MaterialData
{
public:
    std::unique_ptr<int[]> mat_indices;
    std::unique_ptr<Real[]> speculers;
    std::unique_ptr<Real[]> diffuses;
    std::unique_ptr<Real[]> transmits;
    std::unique_ptr<Real[]> emissions;
    std::unique_ptr<Real[]> roughnesss;
    std::unique_ptr<Real[]> nis;
};

template<class Real>
class VertexData
{
public:
    //uvs
    std::unique_ptr<int[]> uv_indices;
    std::unique_ptr<Real[]> uvs;
    //normals
    std::unique_ptr<int[]> normal_indices;
    std::unique_ptr<Real[]> normals;
};

float gamma(float x)
{
    return std::pow(x, 1/2.2f);
}

//Utility Functions!!
void SaveImage(const float *rgb, int width, int height)  // To PNG
{
   std::unique_ptr<unsigned char[]> pixel_colors = std::make_unique<unsigned char[]>(3 * width * height);
  for (int y = 0; y < height; y++) 
  {
    for (int x = 0; x < width; x++) 
    {
      int index = y * width + x;
      pixel_colors[index * 3 + 0] = static_cast<unsigned char>(std::max(0.0f, std::min(gamma(rgb[index * 3 + 0]) * 255.0f, 255.0f)));
      pixel_colors[index * 3 + 1] = static_cast<unsigned char>(std::max(0.0f, std::min(gamma(rgb[index * 3 + 1]) * 255.0f, 255.0f)));
      pixel_colors[index * 3 + 2] = static_cast<unsigned char>(std::max(0.0f, std::min(gamma(rgb[index * 3 + 2]) * 255.0f, 255.0f)));
    }
  }
  stbi_write_png("output.png", width, height, 3, &(pixel_colors[0]), width * 3);
}

class OBJloader
{
public:
    OBJloader(const std::string& filename, const std::string& _mtldir):inputfile(filename), mtldir(_mtldir){}
    std::string inputfile;
    std::string mtldir;
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string warn;
    std::string err;
};

void LoadObj_Single_Object(OBJloader& loader, std::vector<float>& vertices, std::vector<int>& indices, 
                          MaterialData<float>& mat_infos, VertexData<float>& vertex_infos, float scale)
{

    bool ret = tinyobj::LoadObj(&loader.attrib, &loader.shapes, &loader.materials, 
                                &loader.warn, &loader.err, loader.inputfile.c_str(), loader.mtldir.c_str());

    // error handling
    if (!loader.warn.empty()) {
        std::cout << loader.warn << std::endl;
    }

    if (!loader.err.empty()) {
        std::cerr << loader.err << std::endl;
    }

    if (!ret) {
        exit(1);
    }




    // load vertices and vertex info
    vertices.resize(loader.attrib.vertices.size());
    for(size_t i = 0; i < loader.attrib.vertices.size(); ++i)
    {
        vertices[i] = loader.attrib.vertices[i];
        
        if(i % 3 == 0) vertices[i] *= -1.0f;
    }

    // if(loader.attrib.normals.size() > 0)
    // {
    //     vertex_infos.normals.reset(new float[loader.attrib.normals.size()]);
    //     for(int i = 0; i < loader.attrib.normals.size(); ++i)
    //     {
    //         vertex_infos.normals[i] = loader.attrib.normals[i];
    //     }
    // }

    // if(loader.attrib.texcoords.size() > 0)
    // {
    //     vertex_infos.uvs.reset(new float[loader.attrib.texcoords.size()]);
    //     for(int i = 0; i < loader.attrib.texcoords.size(); ++i)
    //     {
    //         vertex_infos.uvs[i] = loader.attrib.texcoords[i];
    //     }
    // }

    int num_faces = 0;
    for(size_t i  = 0; i < loader.shapes.size(); ++i)
    {
        num_faces += static_cast<int>(loader.shapes[i].mesh.indices.size()/3);
    }




    //index setting
    // vertex_infos.normal_indices.reset(new int[num_faces * 3 * 3]);
    // vertex_infos.uv_indices.reset(new int[num_faces * 3 * 2]);
    indices.resize(num_faces * 3);
    mat_infos.mat_indices.reset(new int[num_faces]);
    size_t offset = 0;
    for(size_t i = 0; i < loader.shapes.size(); ++i)
    {
        for(size_t f = 0; f < loader.shapes[i].mesh.indices.size()/3; ++f)
        {
            // vertex_infos.normal_indices[3 * (offset + f) + 0] = 
            //     loader.shapes[i].mesh.indices[3 * f + 0].normal_index;
            // vertex_infos.normal_indices[3 * (offset + f) + 1] = 
            //     loader.shapes[i].mesh.indices[3 * f + 1].normal_index;
            // vertex_infos.normal_indices[3 * (offset + f) + 2] = 
            //     loader.shapes[i].mesh.indices[3 * f + 2].normal_index;

            // vertex_infos.uv_indices[3 * (offset + f) + 0] =
            //     loader.shapes[i].mesh.indices[3 * f + 0].texcoord_index;
            // vertex_infos.uv_indices[3 * (offset + f) + 1] =
            //     loader.shapes[i].mesh.indices[3 * f + 1].texcoord_index;
            // vertex_infos.uv_indices[3 * (offset + f) + 2] =
            //     loader.shapes[i].mesh.indices[3 * f + 2].texcoord_index;

            indices[3 * (offset + f) + 0] = 
                loader.shapes[i].mesh.indices[3 * f + 0].vertex_index;
            indices[3 * (offset + f) + 1] = 
                loader.shapes[i].mesh.indices[3 * f + 1].vertex_index;
            indices[3 * (offset + f) + 2] = 
                loader.shapes[i].mesh.indices[3 * f + 2].vertex_index;

            mat_infos.mat_indices[offset + f] = 
                loader.shapes[i].mesh.material_ids[f];
        }
        offset += loader.shapes[i].mesh.indices.size()/3;
    }




    //material setting
    mat_infos.diffuses.reset(new float[3 * loader.materials.size()]);
    mat_infos.emissions.reset(new float[3 * loader.materials.size()]);
    mat_infos.speculers.reset(new float[3 * loader.materials.size()]);
    mat_infos.transmits.reset(new float[3 * loader.materials.size()]);
    //mat_infos.roughnesss.reset(new float[loader.materials.size()]);
    mat_infos.nis.reset(new float[loader.materials.size()]);
    for(int i = 0; i < loader.materials.size(); ++i)
    {
        //diffuse
        mat_infos.diffuses[3 * i + 0] = loader.materials[i].diffuse[0];
        mat_infos.diffuses[3 * i + 1] = loader.materials[i].diffuse[1];
        mat_infos.diffuses[3 * i + 2] = loader.materials[i].diffuse[2];

        //emission
        mat_infos.emissions[3 * i + 0] = loader.materials[i].emission[0];
        mat_infos.emissions[3 * i + 1] = loader.materials[i].emission[1];
        mat_infos.emissions[3 * i + 2] = loader.materials[i].emission[2];

        //speculer
        mat_infos.speculers[3 * i + 0] = loader.materials[i].specular[0];
        mat_infos.speculers[3 * i + 1] = loader.materials[i].specular[1];
        mat_infos.speculers[3 * i + 2] = loader.materials[i].specular[2];

        //transmit
        mat_infos.transmits[3 * i + 0] = loader.materials[i].transmittance[0];
        mat_infos.transmits[3 * i + 1] = loader.materials[i].transmittance[1];
        mat_infos.transmits[3 * i + 2] = loader.materials[i].transmittance[2];

        //ior
        mat_infos.nis[i] = loader.materials[i].ior;
    }
}

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



Vec3<float> Trace(const float* firstRay_dir, const float* firstRay_origin, const MaterialData<float>& mat_infos, const IBL& ibl,
                  const VertexData<float>& vertex_infos ,const BinaryBVH<float, TriangleData<float>>& bvh, RandomManager& rnd_manager)
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
                Vec3<float> shadoworigin  = hitPos + 0.001f * orienting_normal;
                float shadoworigin_ary[3] = {shadoworigin[0], shadoworigin[1], shadoworigin[2]};
                Hit<float> shadow_hit;
                if(bvh.Traverse(shadow_hit, shadoworigin_ary, shadowdir_ary, 0))
                {
                    //nothing
                }
                else
                {
                    float costerm = orienting_normal.dot(shadowdir);
                    I[0] = I[0] + throughput[0] * mat_infos.diffuses[3 * mat_infos.mat_indices[hit.GetID()] + 0] * Le[0]
                                                * costerm
                                                / (nee_pdf * M_PI);
                    I[1] = I[1] + throughput[1] * mat_infos.diffuses[3 * mat_infos.mat_indices[hit.GetID()] + 1] * Le[1]
                                                * costerm
                                                / (nee_pdf * M_PI);
                    I[2] = I[2] + throughput[2] * mat_infos.diffuses[3 * mat_infos.mat_indices[hit.GetID()] + 2] * Le[2]
                                                * costerm
                                                / (nee_pdf * M_PI);
                }
            }
            #endif

            //BRDF sampling
            float u1 = rnd_manager.GetRND();
            float u2 = rnd_manager.GetRND();
            Vec3<float> wi = randomCosineHemisphere(u1, u2, orienting_normal);

            comingRay = Ray(hitPos + 0.001f * orienting_normal, wi);

            // store throughput
            throughput[0] *= mat_infos.diffuses[3 * mat_infos.mat_indices[hit.GetID()] + 0];
            throughput[1] *= mat_infos.diffuses[3 * mat_infos.mat_indices[hit.GetID()] + 1];
            throughput[2] *= mat_infos.diffuses[3 * mat_infos.mat_indices[hit.GetID()] + 2];

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






int main()
{
    //IBL
    IBL ibl("../../../map_textures/railway_bridges_16k.hdr");

    //shape
    std::vector<float> vertices;
    std::vector<int> indices;
    MaterialData<float> mat_infos;
    VertexData<float> vertex_infos;
    OBJloader load("../../../objects/debug.obj", "../../../objects");
    float scale = 0.5;
    LoadObj_Single_Object(load, vertices, indices, mat_infos, vertex_infos,scale);



    //BVH
    TriangleData<float> sankaku_data(vertices.data(), indices.data());
    BinaryBVH<float, TriangleData<float>> bvh(sankaku_data, indices.size()/3);
    bvh.BuildBVH(Evaluator::SAH);



    //Camera
    float theta = 70.0f * M_PI/180.0f;
    float phi =   135.0f * M_PI/180.0f;
    float r = 9.0f;
    float x = r * std::sin(theta) * std::cos(phi);
    float y = r * std::cos(theta);
    float z = r * std::sin(theta) * std::sin(phi);

    float cameraPos[3] = {x, y, z};
    float cameraForward[3] = {-x/r, -y/r, -z/r};
    PinholeCamera<float> pincam(cameraPos, cameraForward);



    //Screen, Image
    int width = 512, height = 512;
    float* RGB = new float[width * height * 3];
    for(int i = 0; i < width*height*3; i++)
    {
        RGB[i] = 0.0f;
    }
    float screen_height = 2.0f;
    float pixel_size = screen_height/height;


    int samples = 100;

    //Rendering
    std::function<void(const int*, const int*, RandomManager&)> render = 
        [&](const int* upper_left, const int* bottom_right, RandomManager& rnd_manager)
        {
            for(int x = upper_left[0]; x <  bottom_right[0]; ++x)
            {
                for(int y = upper_left[1]; y < bottom_right[1]; ++y)
                {
                    for(int i = 0; i < samples; ++i)
                    {
                        float u = x * pixel_size - 0.5f * pixel_size * width + pixel_size * rnd_manager.GetRND();
                        float v = y * pixel_size - 0.5f * pixel_size * height + pixel_size * rnd_manager.GetRND();
                        float ray_dir[3], ray_origin[3];
                        pincam.CreateFirstRay(u, v, ray_origin, ray_dir);
                        Vec3<float> result = Trace(ray_dir, ray_origin, mat_infos, ibl, vertex_infos, bvh, rnd_manager);
                        RGB[3*(width * y + x) + 0] += result[0]/samples;
                        RGB[3*(width * y + x) + 1] += result[1]/samples;
                        RGB[3*(width * y + x) + 2] += result[2]/samples;
                    }
                }
            }
        };
    ParallelRender paral(render);
    paral.Execute(width, height, 20);
    SaveImage(RGB, width, height);
}