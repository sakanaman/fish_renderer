#include <camera.hpp>
#include <Accelerator.hpp>
#include <tiny_obj_loader.h>
#include <stb_image_write.h>
#include <concurrent.hpp>
#include <Material.hpp>


//Utility Functions!!
void SaveImage(const float *rgb, int width, int height)  // To PNG
{
   std::unique_ptr<unsigned char[]> pixel_colors = std::make_unique<unsigned char[]>(3 * width * height);
  for (int y = 0; y < height; y++) 
  {
    for (int x = 0; x < width; x++) 
    {
      int index = y * width + x;
      pixel_colors[index * 3 + 0] = static_cast<unsigned char>(std::max(0.0f, std::min(rgb[index * 3 + 0] * 255.0f, 255.0f)));
      pixel_colors[index * 3 + 1] = static_cast<unsigned char>(std::max(0.0f, std::min(rgb[index * 3 + 1] * 255.0f, 255.0f)));
      pixel_colors[index * 3 + 2] = static_cast<unsigned char>(std::max(0.0f, std::min(rgb[index * 3 + 2] * 255.0f, 255.0f)));
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

void LoadObj_Single_Object(OBJloader& loader, std::vector<float>& vertices, std::vector<int>& indices, std::vector<MaterialData<float>>& mats, float scale)
{

    bool ret = tinyobj::LoadObj(&loader.attrib, &loader.shapes, &loader.materials, 
                                &loader.warn, &loader.err, loader.inputfile.c_str(), loader.mtldir.c_str());

    if (!loader.warn.empty()) {
        std::cout << loader.warn << std::endl;
    }

    if (!loader.err.empty()) {
        std::cerr << loader.err << std::endl;
    }

    if (!ret) {
        exit(1);
    }

    vertices.resize(loader.attrib.vertices.size());
    for(size_t i = 0; i < loader.attrib.vertices.size(); ++i)
    {
        vertices[i] = loader.attrib.vertices[i];
        //if(i % 3 == 2) vertices[i] *= -1.0f;
    }

    for(size_t s = 0; s < loader.shapes.size(); ++s)
    {
        size_t index_offset = 0;
        std::cout << loader.shapes[s].mesh.num_face_vertices.size() << std::endl;
        for(size_t f = 0; f < loader.shapes[s].mesh.num_face_vertices.size(); ++f)
        {
            int fv = loader.shapes[s].mesh.num_face_vertices[f];
            for(size_t v = 0; v < fv; ++v)
            {
                tinyobj::index_t idx = loader.shapes[s].mesh.indices[index_offset + v];
                indices.push_back(idx.vertex_index);
            }
            index_offset += fv;
        }
    }
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



Vec3<float> Trace(const float* firstRay_dir, const float* firstRay_origin, const std::vector<MaterialData<float>>& mats,
                  const BinaryBVH<float, TriangleData<float>>& bvh, pcg32_random_t* rng)
{
    Vec3<float> throughput(1.0f, 1.0f, 1.0f);
    Vec3<float> I(0.0f, 0.0f, 0.0f);
    Ray comingRay(firstRay_origin, firstRay_dir);
}






int main()
{
    //shape
    std::vector<float> vertices;
    std::vector<int> indices;
    std::vector<MaterialData<float>> mats;
    OBJloader load("../../objects/mokey_bunny_plane.obj", "../../objects");
    float scale = 0.5;
    LoadObj_Single_Object(load, vertices, indices, mats,scale);



    //BVH
    TriangleData<float> sankaku_data(vertices.data(), indices.data());
    BinaryBVH<float, TriangleData<float>> bvh(sankaku_data, indices.size()/3);
    bvh.BuildBVH(Evaluator::SAH);



    //Camera
    float cameraPos[3] = {0.0f, 2.0f, 9.0f};
    float cameraForward[3] = {0.0f, 0.0f, -1.0f};
    PinholeCamera<float> pincam(cameraPos, cameraForward);



    //Screen, Image
    int width = 1024, height = 1024;
    float* RGB = new float[width * height * 3];
    for(int i = 0; i < width*height*3; i++)
    {
        RGB[i] = 0.0f;
    }
    float screen_height = 2.0f;
    float pixel_size = screen_height/height;




    //Rendering
    std::function<void(const int*, const int*, pcg32_random_t* rng)> render = 
        [&](const int* upper_left, const int* bottom_right, pcg32_random_t* rng)
        {
            for(int x = upper_left[0]; x <  bottom_right[0]; ++x)
            {
                for(int y = upper_left[1]; y < bottom_right[1]; ++y)
                {
                    float u = x * pixel_size - 0.5f * pixel_size * width;
                    float v = y * pixel_size - 0.5f * pixel_size * height;
                    float ray_dir[3], ray_origin[3];
                    pincam.CreateFirstRay(u, v, ray_origin, ray_dir);
                    Vec3<float> result = Trace(ray_dir, ray_origin, mats, bvh, rng);
                    RGB[3*(width * y + x) + 0] = result[0];
                    RGB[3*(width * y + x) + 1] = result[1];
                    RGB[3*(width * y + x) + 2] = result[2];
                }
            }
        };
    ParallelRender paral(render);
    paral.Execute(width, height, 8);
    SaveImage(RGB, width, height);
}