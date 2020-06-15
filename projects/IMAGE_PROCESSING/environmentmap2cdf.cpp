#include <stb_image.h>
#include <stb_image_write.h>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <memory>
#include <random.hpp>

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

int main()
{
    std::string filename = "../../../map_textures/Tokyo_BigSight_3k.hdr";
    int width, height, channel;
    float* data = stbi_loadf(filename.c_str(), &width, &height, &channel, 0);



    // make 1D map
    int width1 = width, height1 = height;
    float* RGB1 = new float[width1 * height1 * 3];
    for(int x = 0; x < width1; ++x)
    {
        for(int y = 0; y < height1; ++y)
        {
            float r = data[3 * (width * y + x) + 0];
            float g = data[3 * (width * y + x) + 1];
            float b = data[3 * (width * y + x) + 2];
            
            float L = 0.2126 * r + 0.7152 * g + 0.0722 * b;
            RGB1[3 * (width * y + x) + 0] = L;
            RGB1[3 * (width * y + x) + 1] = L;
            RGB1[3 * (width * y + x) + 2] = L;
        }
    }
    stbi_write_hdr("environmentL_1D.hdr", width1, height1, 3, RGB1);




    //Make PDF map
    float angles[height];
    for(int i = 0; i < height; i++)
    {
        angles[i] = (2.0f * i + 1) * M_PI / (2.0f * height);
    }

    float sinthetas[height];
    for(int i = 0; i < height; i++)
    {
        sinthetas[i] = std::sin(angles[i]);
    }

    std::vector<std::vector<float>> PVDists(width, std::vector<float>(height));
    std::vector<float> PUDists(width);

    for(int x = 0; x < width; ++x)
    {
        //y = 0
        float r = data[3 * (width * 0 + x) + 0];
        float g = data[3 * (width * 0 + x) + 1];
        float b = data[3 * (width * 0 + x) + 2];
        float L = 0.2126 * r + 0.7152 * g + 0.0722 * b;
        PVDists.at(x).at(0) = sinthetas[0] * L;

        //y = 1, 2,...,height-1
        for(int y = 1; y < height; ++y)
        {
            r = data[3 * (width * y + x) + 0];
            g = data[3 * (width * y + x) + 1];
            b = data[3 * (width * y + x) + 2];
            L = 0.2126 * r + 0.7152 * g + 0.0722 * b;

            PVDists.at(x).at(y) = PVDists.at(x).at(y-1) + L * sinthetas[y];
        }

        if(x == 0)
        {
            PUDists.at(x) = PVDists.at(x).at(height - 1);
        }
        else
        {
            PUDists.at(x) = PUDists.at(x - 1) + PVDists.at(x).at(height - 1);
        }
    }

    float dt = M_PI/height;
    float invPdfNorm = 2.0f * M_PI * M_PI / (width * height);
    
    float* RGB2 = new float[width * height * 3];
    for(int u = 0; u < width; ++u)
    {
        for(int v = 0; v < height; ++v)
        {
            float pdfU = (u == 0) ? PUDists.at(0) : PUDists.at(u) - PUDists.at(u-1);
            pdfU /= PUDists.at(width - 1);
            float pdfV = (v == 0) ? PVDists.at(u).at(0) : PVDists.at(u).at(v) - PVDists.at(u).at(v-1);
            pdfV /= PVDists.at(u).at(height - 1);

            float theta = dt * 0.5 + dt * v;
            float PDF = pdfV * pdfU / (invPdfNorm * std::sin(theta));

            RGB2[3 * (width * v + u) + 0] = PDF;
            RGB2[3 * (width * v + u) + 1] = PDF;
            RGB2[3 * (width * v + u) + 2] = PDF;
        }
    }
    stbi_write_hdr("environmentPDF_1D.hdr", width, height, 3, RGB2);


    //Let's sample!
    int* count_list = new int[width * height];
    for(int i = 0; i < width * height; ++i)
    {
        count_list[i] = 0;
    }

    pcg32_random_t rng;
    int samples = 10000000;
    for(int i = 0; i < samples; ++i)
    {
        int u, v;

        float maxUval = PUDists[width - 1];
        auto pUPos = std::lower_bound(PUDists.begin(), PUDists.end(), maxUval * rng.rnd());
        u = pUPos - PUDists.begin();

        float maxVval = PVDists.at(u).at(height - 1);
        auto pVPos = std::lower_bound(PVDists.at(u).begin(), PVDists.at(u).end(), maxVval * rng.rnd());
        v = pVPos - PVDists.at(u).begin();

        count_list[width * v + u] += 1;
    }

    int serch_most_sampled = 0;
    int most_sample_u = 0;
    for(int i = 0; i < width*height; ++i)
    {
        if(serch_most_sampled < count_list[i])
        {
            serch_most_sampled = count_list[i];
            most_sample_u = i;
        }
    }
    assert(serch_most_sampled > 0);
    std::cout <<"[" << most_sample_u << "]: " << serch_most_sampled << std::endl;

    //input hdr
    float* RGB3 = new float[width * height * 3];
    for(int x = 0; x < width; ++x)
    {
        for(int y = 0; y < height; ++y)
        {
            RGB3[3 * (width * y + x) + 0] = float(count_list[width*y + x])/serch_most_sampled;
            RGB3[3 * (width * y + x) + 1] = float(count_list[width*y + x])/serch_most_sampled;
            RGB3[3 * (width * y + x) + 2] = float(count_list[width*y + x])/serch_most_sampled;
            //std::cout << float(count_list[width*y + x])/serch_most_sampled << std::endl;
        }
    }
    SaveImage(RGB3, width, height);
}