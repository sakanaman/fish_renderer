#ifndef SCENE_HPP
#define SCENE_HPP

#include "api.hpp"
#include <memory>
#include <iostream>

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

template<class Real>
class SceneData
{
public:
private:
    VertexData<Real> v_data;
    MaterialData<Real> mat_data;
};

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


#endif