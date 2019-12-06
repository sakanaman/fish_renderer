#include "IBL.hpp"
#include <iostream>
#include <cmath>

IBL::IBL(std::string _filename):name(_filename){
    std::cout << "Load IBL: " << name.c_str() << std::endl;
    data = stbi_loadf(_filename.c_str(), &w, &h, &channel, 0);
    std::cout << channel << std::endl;
    if(data){   
        std::cout << "load success!!" << std::endl;
    }
    std::cout << "data[0]" << data[13] << std::endl;

    std::cout << "Make IBL sampler" << std::endl;
    PVDists = std::vector<std::vector<float>>(w, std::vector<float>(h));
    PUDists =  std::vector<float>(w);

    float angles[h];
    for(int i = 0; i < h; i++)
    {
        angles[i] = (2.0f * i + 1) * M_PI / (2.0f * h);
    }

    float sinthetas[h];
    for(int i = 0; i < h; i++)
    {
        sinthetas[i] = std::sin(angles[i]);
    }

    for(int x = 0; x < w; ++x)
    {
        //y = 0
        float r = data[3 * (w * 0 + x) + 0];
        float g = data[3 * (w * 0 + x) + 1];
        float b = data[3 * (w * 0 + x) + 2];
        float L = 0.2126 * r + 0.7152 * g + 0.0722 * b;
        PVDists.at(x).at(0) = sinthetas[0] * L;

        //y = 1, 2,...,height-1
        for(int y = 1; y < h; ++y)
        {
            r = data[3 * (w * y + x) + 0];
            g = data[3 * (w * y + x) + 1];
            b = data[3 * (w * y + x) + 2];
            L = 0.2126 * r + 0.7152 * g + 0.0722 * b;

            PVDists.at(x).at(y) = PVDists.at(x).at(y-1) + L * sinthetas[y];
        }

        if(x == 0)
        {
            PUDists.at(x) = PVDists.at(x).at(h - 1);
        }
        else
        {
            PUDists.at(x) = PUDists.at(x - 1) + PVDists.at(x).at(h - 1);
        }
    }
    std::cout << "finish make IBL sampler" << std::endl;
}

void IBL::xy2uv(const int x, const int y, float* u, float* v) const 
{
    *u = float(x)/w;
    *v = float(y)/h;
}

void IBL::uv2xy(int* x, int* y, const float u, const float v) const
{
    *x = static_cast<int>(w * u);
    *y = static_cast<int>(h * v);
}

void IBL::angle2uv(float* u, float* v, const float theta, const float phi) const 
{
    *v = theta / M_PI;
    *u = phi / (2 * M_PI);
}

void IBL::GetLe(float* Le, const float theta, const float phi) const
{
    float modefy_phi = 2.0f*M_PI - phi;//modify inversion
    float u, v;
    angle2uv(&u, &v, theta, modefy_phi);
    int x, y;
    uv2xy(&x, &y, u, v);
    Le[0] = data[3 * (w * y + x) + 0];
    Le[1] = data[3 * (w * y + x) + 1];
    Le[2] = data[3 * (w * y + x) + 2];
}

float IBL::sample(const float r1, const float r2, float* phi, float* theta, float* Le) const//return PDF
{
    int u, v;

    float maxUval = PUDists[w - 1];
    auto pUPos = std::lower_bound(PUDists.begin(), PUDists.end(), maxUval * r1);
    u = pUPos - PUDists.begin();

    float maxVval = PVDists.at(u).at(h - 1);
    auto pVPos = std::lower_bound(PVDists.at(u).begin(), PVDists.at(u).end(), maxVval * r2);
    v = pVPos - PVDists.at(u).begin();

    float dt = M_PI/h;
    float invPdfNorm = 2.0f * M_PI * M_PI / (w * h);
    float pdfU = (u == 0) ? PUDists.at(0) : PUDists.at(u) - PUDists.at(u-1);
    pdfU /= PUDists.at(w - 1);
    float pdfV = (v == 0) ? PVDists.at(u).at(0) : PVDists.at(u).at(v) - PVDists.at(u).at(v-1);
    pdfV /= PVDists.at(u).at(h - 1);

    Le[0] = data[3 * (w * v + u) + 0];
    Le[1] = data[3 * (w * v + u) + 1];
    Le[2] = data[3 * (w * v + u) + 2];
    *theta = dt * 0.5 + dt * v;
    *phi = 2.0f * M_PI - float(u)/w * 2.0f * M_PI;//modify inversion
    float PDF = pdfV * pdfU / (invPdfNorm * std::sin(*theta));
    return PDF;
}