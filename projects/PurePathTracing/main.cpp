#define NEE
#include "raytrace.hpp"

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



int main()
{
    //IBL
    IBL ibl("../../../map_textures/PaperMill_E_3k.hdr");

    //shape
    std::vector<float> vertices;
    std::vector<int> indices;
    SceneData<float> scenedata;
    OBJloader load("../../../objects/debug.obj", "../../../objects");
    float scale = 0.5;
    LoadObj_Single_Object(load, vertices, indices, scenedata,scale);


    //Shader
    Shader<float> shader = DiffuseShader{scenedata};


    //BVH
    TriangleData<float> sankaku_data(vertices.data(), indices.data());
    BinaryBVH<float, TriangleData<float>> bvh(sankaku_data, indices.size()/3);
    bvh.BuildBVH(Evaluator::SAH);


    //Camera
    float theta = 70.0f * M_PI/180.0f;
    float phi =   270.0f * M_PI/180.0f;
    float r = 6.0f;
    float x = r * std::sin(theta) * std::cos(phi);
    float y = r * std::cos(theta);
    float z = r * std::sin(theta) * std::sin(phi);

    float cameraPos[3] = {x, y + 2.0f, z};
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
                        //Vec3<float> result = Trace(ray_dir, ray_origin, scenedata, ibl, bvh, rnd_manager);
                        Vec3<float> result = Trace_Test(ray_dir, ray_origin, shader, ibl, bvh, rnd_manager);
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