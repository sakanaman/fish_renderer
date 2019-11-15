#include "Accelerator.hpp"
#include <cmath>
#include <bitset>
#include <iostream>

//Hit class implementation
template<class Real>
Hit<Real>::Hit()
{
    t = 1000000;
}


template<class Real>
Real Hit<Real>::GetT() const
{
    return t;
}

template<class Real>
void Hit<Real>::SetT(const Real _t)
{
    t = _t;
}

template<class Real>
Real* Hit<Real>::GetPos()
{
    return point;
}

template<class Real>
Real* Hit<Real>::GetNg()
{
    return normal;
}

template<class Real>
int Hit<Real>::GetID() const
{
    return geomID;
}

template<class Real>
void Hit<Real>::SetPos(const Real* pos)
{   
    point[0] = pos[0];
    point[1] = pos[1];
    point[2] = pos[2];
}

template<class Real>
void Hit<Real>::SetNg(const Real* n)
{
    normal[0] = n[0];
    normal[1] = n[1];
    normal[2] = n[2];
}

template<class Real>
void Hit<Real>::SetID(const int num)
{
    geomID = num;
}

template class Hit<float>;
template class Hit<double>;

//utility enum
enum class Axis
{
    X,
    Y,
    Z
};

enum class LR
{
    Left,
    Right
};

//AABB class implementation
template<class Real>
AABB<Real>::AABB(const Real* _max, const Real* _min)
{
    max[0] = _max[0];
    max[1] = _max[1];
    max[2] = _max[2];

    min[0] = _min[0];
    min[1] = _min[1];
    min[2] = _min[2];
}

template<class Real>
AABB<Real>::AABB()
{}

template<class Real>
void AABB<Real>::Setter(const Real* _max, const Real* _min)
{
    max[0] = _max[0];
    max[1] = _max[1];
    max[2] = _max[2];

    min[0] = _min[0];
    min[1] = _min[1];
    min[2] = _min[2];
}

template<class Real>
bool AABB<Real>::intersect(const Real* ray_origin, const Real* ray_dir) const {
    Real t_max = std::numeric_limits<Real>::max();
    Real t_min = std::numeric_limits<Real>::lowest();

    for(int i = 0; i < 3; ++i) {
        Real t1 = (min[i] - ray_origin[i]) / ray_dir[i];
        Real t2 = (max[i] - ray_origin[i]) / ray_dir[i];
        
        Real t_near = std::min(t1, t2);
        Real t_far = std::max(t1, t2);
        
        t_max = std::min(t_max, t_far);
        t_min = std::max(t_min, t_near);
    
        if (t_min > t_max) return false;
    }
    return true;
}

template<class Real>
AABB<Real> AABB<Real>::Union(const AABB& r) const
{
    AABB result;

    for(int i = 0; i < 3; ++i) {
        result.max[i] =  std::max(max[i], r.max[i]);
        result.min[i] =  std::min(min[i], r.min[i]);
    }
    return result;
}

template<class Real>
AABB<Real> AABB<Real>::Union(const Real a[3]) const
{
    AABB result;

    for(int i = 0; i < 3; ++i) {
        result.max[i] = std::max(max[i], a[i]);
        result.min[i] = std::min(min[i], a[i]);
    }
    return result;
}

template<class Real>
Real AABB<Real>::GetMax(const Axis axis) const
{
    return max[static_cast<int>(axis)];
}

template<class Real>
Real AABB<Real>::GetMin(const Axis axis) const
{ 
    return min[static_cast<int>(axis)];
}

template<class Real>
Real AABB<Real>::AreaAABB() const
{
    Real result;
    result = 2*((max[0] - min[0])*(max[1] - min[1]))
            +2*((max[1] - min[1])*(max[2] - min[2]))
            +2*((max[2] - min[2])*(max[0] - min[0]));
    return result;
}

template class AABB<float>;
template class AABB<double>;

//BinaryBVH`s Node implementation

template<class Real>
void BBVHNode<Real>::SetSplitDimention(const Axis axis)
{
    switch (axis)
    {
        case Axis::X :
        {
            // nothing
            break;
        }
        case Axis::Y :
        {
            childinfo_left |= (1 << 29);
            break;   
        }
        case Axis::Z :
        {
            childinfo_left |= (1 << 30);
            break;
        }
    }
    childinfo_left |= (1 << 31);
}

template<class Real>
void BBVHNode<Real>::SetChildIndex(const LR lr, const int index)
{
    switch (lr)
    {
        case LR::Left :
        {
            childinfo_left |= index;
            break;
        }
        case LR::Right :
        {
            childinfo_right |= index;
            break;
        }
    }
}

template<class Real>
void BBVHNode<Real>::SetLeafNode(const int* range)
{
    childinfo_left = range[0];
    childinfo_right = range[1];
}

template<class Real>
void BBVHNode<Real>::SetAABB(const AABB<Real>& _aabb)
{
    aabb = _aabb;
}

template<class Real>
AABB<Real> BBVHNode<Real>::GetAABB() const
{
    return aabb;
}

template<class Real>
int BBVHNode<Real>::GetBegin() const
{
    return childinfo_left;
}

template<class Real>
int BBVHNode<Real>::GetEnd() const
{
    return childinfo_right;
}

template<class Real>
bool BBVHNode<Real>::isLeaf() const
{
    return !(childinfo_left & (1 << 31));
}

template<class Real>
int BBVHNode<Real>::GetChild(const LR lr)
{
    switch(lr)
    {
        case LR::Left :
        {
            return ((childinfo_left << 3) >> 3);
            break;
        }
        case LR::Right :
        {
            return ((childinfo_right << 3) >> 3);
            break;
        }
    }
}

template class BBVHNode<float>;
template class BBVHNode<double>;


//BinaryBVH implementation
template<class Real, class ShapeData>
BinaryBVH<Real, ShapeData>::BinaryBVH(const ShapeData& _shapedata, const int face_num)
{
    faces = std::vector<int>(face_num);
    for(int i = 0; i < face_num; ++i)
    {
        // set geomID to faces
        faces[i] = i;
    }
    shapedata = _shapedata;
    aabbs = std::make_unique<AABB<Real>[]>(face_num);
    for(int i = 0; i < face_num; ++i)
    {
        shapedata.makeAABB(i, aabbs[i]);
    }
    bvh_nodes = std::make_unique<BBVHNode<Real>[]>(3 * face_num);
}

template<class Real>
class BucketInfo {
public:
    int count = 0;
    AABB<Real> aabb;
};

template class BucketInfo<float>;
template class BucketInfo<double>;

template<class Real, class ShapeData>
int BinaryBVH<Real, ShapeData>::PartitionSAH(const int* range, const AABB<Real>& bigaabb)
{
    //std::cout << "this task:  [" << range[0] << ", " << range[1] << ")" << std::endl;
    if(range[1] - range[0] == 1)
    {
        bvh_nodes[range[2]].SetLeafNode(range);
        return -1;
    }
    else
    {
        AABB<Real> centoroidAABB;
        for(int i = range[0]; i < range[1]; ++i)
        {
            Real centroid[3];
            for(int j = 0; j < 3; ++j)
            {
                centroid[j] =( aabbs[faces[i]].GetMax(static_cast<Axis>(j))
                              +aabbs[faces[i]].GetMin(static_cast<Axis>(j))
                             ) * 0.5;
            }
            centoroidAABB = centoroidAABB.Union(centroid);
        }

        //choose dimention
        int dim;
        Real s = -1.0;
        for(int i = 0; i < 3; ++i)
        {
            Real max_num = centoroidAABB.GetMax(static_cast<Axis>(i));
            Real min_num = centoroidAABB.GetMin(static_cast<Axis>(i));
            if(s < max_num - min_num)
            {
                s = max_num - min_num;
                dim = i; 
            }
        }

        if(centoroidAABB.GetMax(static_cast<Axis>(dim)) == centoroidAABB.GetMin(static_cast<Axis>(dim)))
        {
            bvh_nodes[range[2]].SetLeafNode(range);
            return -1;
        }
        else
        {
            if(range[1] - range[0] <= 4)
            {
                int mid = (range[1] + range[0])/2;
                std::nth_element(faces.begin() + range[0], faces.begin() + mid, faces.begin() + range[1],
                                 [&](const int a, const int b) {
                                    return aabbs[a].GetMax(static_cast<Axis>(dim)) + aabbs[a].GetMin(static_cast<Axis>(dim))
                                         < aabbs[b].GetMax(static_cast<Axis>(dim)) + aabbs[b].GetMin(static_cast<Axis>(dim));   
                                 });
                bvh_nodes[range[2]].SetSplitDimention(static_cast<Axis>(dim));
                return mid;
            }
            else
            {
                const int nBuckets = 12;
                BucketInfo<Real> buckets[nBuckets];
                for(int i = range[0]; i < range[1]; ++i)
                {
                    int b = nBuckets * 
                            ((aabbs[faces[i]].GetMax(static_cast<Axis>(dim)) + aabbs[faces[i]].GetMin(static_cast<Axis>(dim)))*0.5 - centoroidAABB.GetMin(static_cast<Axis>(dim)))
                            / (centoroidAABB.GetMax(static_cast<Axis>(dim)) - centoroidAABB.GetMin(static_cast<Axis>(dim)));
                    if(b == nBuckets) b = nBuckets - 1;
                    buckets[b].count++;
                    buckets[b].aabb = (buckets[b].aabb).Union(aabbs[faces[i]]);
                }

                Real cost[nBuckets -1];
                for(int i = 0; i < nBuckets -1; ++i)
                {
                    AABB<Real> b0, b1;
                    int count0 = 0, count1 = 0;
                    for(int j = 0; j <= i; j++)
                    {
                        b0 = b0.Union(buckets[j].aabb);
                        count0 += buckets[j].count;
                    }
                    for(int j = i+1; j < nBuckets; j++)
                    {
                        b1 = b1.Union(buckets[j].aabb);
                        count1 += buckets[j].count;
                    }
                    cost[i]= 0.125 + (count0 * b0.AreaAABB() + 
                                      count1 * b1.AreaAABB()) / bigaabb.AreaAABB();
                    assert(!std::isnan(cost[i]));
                }

                Real minCost = cost[0];
                int minCostSplitBucket = 0;
                for(int i = 0; i < nBuckets - 1; ++i)
                {
                    if(cost[i] < minCost)
                    {
                        minCost = cost[i];
                        minCostSplitBucket = i;
                    }
                }

                Real leafCost = range[1] - range[0];
                if(minCost < leafCost)
                {
                    bvh_nodes[range[2]].SetSplitDimention(static_cast<Axis>(dim));
                    auto pos = std::partition(faces.begin() + range[0], faces.begin() + range[1], 
                                            [&](const int geomID)
                                            {
                                                int b = nBuckets * 
                                                        ((aabbs[geomID].GetMax(static_cast<Axis>(dim)) + aabbs[geomID].GetMin(static_cast<Axis>(dim))) * 0.5 - centoroidAABB.GetMin(static_cast<Axis>(dim)))
                                                        / (centoroidAABB.GetMax(static_cast<Axis>(dim)) - centoroidAABB.GetMin(static_cast<Axis>(dim)));
                                                if(b == nBuckets) b = nBuckets - 1;
                                                return b <= minCostSplitBucket;
                                            });
                    return pos - faces.begin();
                }
                else
                {
                    bvh_nodes[range[2]].SetLeafNode(range);
                    return -1;
                }
                
            }
            
        }
        
    }
};

template<class Real, class ShapeData>
void BinaryBVH<Real, ShapeData>::BuildBVH(const Evaluator eval)
{
    int nodecount = 0;
    int Range[1024*3];
    int remaintasks;
    //initialize task
    remaintasks = 1;
    Range[0] = 0;
    Range[1] = faces.size();
    Range[2] = 0;

    int mytask[3];

    while(1) 
    {
        {//fetch task
            if(remaintasks == 0)
            {
                break;
            }
            --remaintasks;
            mytask[0] = Range[remaintasks*3 + 0];
            mytask[1] = Range[remaintasks*3 + 1];
            mytask[2] = Range[remaintasks*3 + 2];
        }

        for(int i = mytask[0]; i < mytask[1]; ++i)
        {
            bvh_nodes[mytask[2]].SetAABB((bvh_nodes[mytask[2]].GetAABB()).Union(aabbs[faces[i]]));
        }

        AABB<Real> bigaabb = bvh_nodes[mytask[2]].GetAABB();
        int bestsplit = PartitionSAH(mytask, bigaabb);
        {
            if(bestsplit != -1)
            {
                //std::cout << bestsplit << std::endl;
                Range[remaintasks*3 + 0] = mytask[0];
                Range[remaintasks*3 + 1] = bestsplit;
                Range[remaintasks*3 + 2] = nodecount + 1;
                bvh_nodes[mytask[2]].SetChildIndex(LR::Left, nodecount + 1);

                Range[remaintasks*3 + 3] = bestsplit;
                Range[remaintasks*3 + 4] = mytask[1];
                Range[remaintasks*3 + 5] = nodecount + 2;
                bvh_nodes[mytask[2]].SetChildIndex(LR::Right, nodecount + 2);
                remaintasks += 2;
                nodecount += 2;
            }
            else
            {
                //std::cout << "non split" << std::endl;
            }
        }
    }
    std::cout << "remmainthreads: " << remaintasks << std::endl;
    std::cout << "nodecount: " << nodecount << std::endl;
    std::cout << "root node's left child: " << std::bitset<32>(bvh_nodes[0].GetChild(LR::Left))<< std::endl;
    std::cout << "root node's right child: " << std::bitset<32>(bvh_nodes[0].GetChild(LR::Right))<< std::endl;
}

template<class Real, class ShapeData>
bool BinaryBVH<Real, ShapeData>::Traverse(Hit<Real>& hit, const Real* ray_origin, const Real* ray_dir, const int index) const
{
    bool is_hit_aabb = (bvh_nodes[index].GetAABB()).intersect(ray_origin, ray_dir);
    if(!is_hit_aabb){
        return false;
    }
    else
    {
        if(bvh_nodes[index].isLeaf())
        {
            Hit<Real> hit_each;
            bool is_hit = false;
            for(int i = bvh_nodes[index].GetBegin(); i < bvh_nodes[index].GetEnd(); ++i)
            {
                if(shapedata.intersect(faces[i], ray_origin, ray_dir, hit_each))
                {
                    if(hit_each.GetT() < hit.GetT())
                    {
                        hit = hit_each;
                    }
                    is_hit = true;
                }
            }
            return is_hit;
        }
        else
        {
            bool is_hit1 = Traverse(hit, ray_origin, ray_dir, bvh_nodes[index].GetChild(LR::Left));
            bool is_hit2 = Traverse(hit, ray_origin, ray_dir, bvh_nodes[index].GetChild(LR::Right));
            return is_hit1 || is_hit2;
        }
        
    }

    // Hit<Real> hit_each;
    // bool is_hit = false;
    // for(int i = 0; i < faces.size(); ++i)
    // {
    //     if(shapedata.intersect(faces[i], ray_origin, ray_dir, hit_each))
    //     {
    //         if(hit_each.GetT() < hit.GetT())
    //         {
    //             hit = hit_each;
    //         }
    //         is_hit = true;
    //     }
    // }
    // return is_hit;
    
}


//shapedata implementation
template<class Real>
SphereData<Real>::SphereData(Real* rad_cents):rad_cents(rad_cents){}

template<class Real>
SphereData<Real>::SphereData(){}

template<class Real>
bool SphereData<Real>::intersect(const int geomID, const Real* ray_origin, const Real* ray_dir, Hit<Real>& hit) const//reference from smallpt
{
    Real radius   = rad_cents[4 * geomID + 0];
    Real center_x = rad_cents[4 * geomID + 1];
    Real center_y = rad_cents[4 * geomID + 2];
    Real center_z = rad_cents[4 * geomID + 3];

    Real d_norm = std::sqrt(ray_dir[0]*ray_dir[0] 
                          + ray_dir[1]*ray_dir[1] 
                          + ray_dir[2]*ray_dir[2]);
    Real oc_norm = std::sqrt( (ray_origin[0] - center_x) * (ray_origin[0] - center_x) 
                            + (ray_origin[1] - center_y) * (ray_origin[1] - center_y)
                            + (ray_origin[2] - center_z) * (ray_origin[2] - center_z));
    Real a = d_norm * d_norm;
    Real b = 2 * ( ray_dir[0] * (ray_origin[0] - center_x)
                 + ray_dir[1] * (ray_origin[1] - center_y)
                 + ray_dir[2] * (ray_origin[2] - center_z));
    Real c =  oc_norm * oc_norm - radius * radius;
    Real d = b*b - 4*a*c;

    if (d < 0)
    {
        return false;
    }
     

    Real t1 = (-b - sqrt(d)) / (2 * a);
    Real t2 = (-b + sqrt(d)) / (2 * a);

    Real t = t1;
    if (t <= 1e-6)
    {
        t = t2;
        if (t <= 1e-6)
        {
            return false;
        }
    }
    hit.SetT(t);
    Real hitPos[3] = {ray_origin[0] + t * ray_dir[0], 
                      ray_origin[1] + t * ray_dir[1], 
                      ray_origin[2] + t * ray_dir[2]};
    hit.SetPos(hitPos);
    Real hitNormal[3] = {hitPos[0] - center_x, 
                         hitPos[1] - center_y, 
                         hitPos[2] - center_z};
    Real normal_len = std::sqrt( hitNormal[0] * hitNormal[0]+
                                 hitNormal[1] * hitNormal[1]+
                                 hitNormal[2] * hitNormal[2]);
    hitNormal[0] /= normal_len;
    hitNormal[1] /= normal_len;
    hitNormal[2] /= normal_len;
    hit.SetNg(hitNormal);
    hit.SetID(geomID);
    return true;
}

template<class Real>
void SphereData<Real>::makeAABB(const int geomID, AABB<Real>& aabb) const
{
    Real radius   = rad_cents[4 * geomID + 0];
    Real center_x = rad_cents[4 * geomID + 1];
    Real center_y = rad_cents[4 * geomID + 2];
    Real center_z = rad_cents[4 * geomID + 3];

    Real max[] = {center_x + radius, center_y + radius, center_z + radius};
    Real min[] = {center_x - radius, center_y - radius, center_z - radius};
    aabb.Setter(max, min);
}

template<class Real>
TriangleData<Real>::TriangleData(Real* _vertices, int* _indices)
                                :vertices(_vertices), indices(_indices){}

template<class Real>
TriangleData<Real>::TriangleData(){};


template<class Real>
bool TriangleData<Real>::intersect(const int geomID, const Real* ray_origin , const Real* ray_dir, Hit<Real>& hit) const
{
    Real v0[3], v1[3], v2[3];

    int index[3] = {indices[3 * geomID + 0], 
                    indices[3 * geomID + 1],
                    indices[3 * geomID + 2]};

    v0[0] = vertices[3 * index[0] + 0];
    v0[1] = vertices[3 * index[0] + 1];
    v0[2] = vertices[3 * index[0] + 2];

    v1[0] = vertices[3 * index[1] + 0];
    v1[1] = vertices[3 * index[1] + 1];
    v1[2] = vertices[3 * index[1] + 2];

    v2[0] = vertices[3 * index[2] + 0];
    v2[1] = vertices[3 * index[2] + 1];
    v2[2] = vertices[3 * index[2] + 2];

    Real e1[3] = {v1[0] - v0[0],
                   v1[1] - v0[1],
                   v1[2] - v0[2]};
    Real e2[3] = {v2[0] - v0[0],
                   v2[1] - v0[1],
                   v2[2] - v0[2]};
    Real d[3] = {ray_dir[0], 
                  ray_dir[1], 
                  ray_dir[2]};
    
    Real r[3] = {ray_origin[0] - v0[0],
                  ray_origin[1] - v0[1],
                  ray_origin[2] - v0[2]};

    Real alpha[3] = {d[1] * e2[2] - d[2] * e2[1],
                      d[2] * e2[0] - d[0] * e2[2],
                      d[0] * e2[1] - d[1] * e2[0]};
    Real beta[3] =  {r[1] * e1[2] - r[2] * e1[1],
                      r[2] * e1[0] - r[0] * e1[2],
                      r[0] * e1[1] - r[1] * e1[0]};
    Real denominator = alpha[0] * e1[0] + 
                        alpha[1] * e1[1] + 
                        alpha[2] * e1[2];

    if(std::abs(denominator) < 1e-7)
    {
        return false;
    }
    else
    {
        Real _t = 1.0f/denominator * (beta[0] * e2[0] + 
                                       beta[1] * e2[1] + 
                                       beta[2] * e2[2]);
        Real x = 1.0f/denominator * (alpha[0] * r[0] + 
                                      alpha[1] * r[1] + 
                                      alpha[2] * r[2]);
        Real y = 1.0f/denominator * (beta[0] * d[0] + 
                                      beta[1] * d[1] + 
                                      beta[2] * d[2]);

        if(_t > 1e-6 && 0 <= x && x <= 1 && 0 <= y && y <= 1 && 0 <= x + y && x + y <= 1)
        {
            hit.SetT(_t);
            hit.SetID(geomID);

            Real hitNormal[3] = {e2[1]*e1[2] - e2[2]*e1[1],
                                  e2[2]*e1[0] - e2[0]*e1[2],
                                  e2[0]*e1[1] - e2[1]*e1[0]};
            Real ng_len = std::sqrt(hitNormal[0] * hitNormal[0] + 
                                     hitNormal[1] * hitNormal[1] +
                                     hitNormal[2] * hitNormal[2]);
            hitNormal[0] /= -1.0 * ng_len;
            hitNormal[1] /= -1.0 * ng_len;
            hitNormal[2] /= -1.0 * ng_len;
            hit.SetNg(hitNormal);

            Real hitPos[3] = {ray_origin[0] + _t * ray_dir[0],
                               ray_origin[1] + _t * ray_dir[1],
                               ray_origin[2] + _t * ray_dir[2]};
            hit.SetPos(hitPos);
            
            return true;
        }
        else
        {
            return false;
        }
    }
}

template<class Real>
void TriangleData<Real>::makeAABB(const int geomID, AABB<Real>& aabb) const
{
    Real v0[3], v1[3], v2[3];
    //std::cout << indices[3 * geomID + 0] << std::endl;
    int index[3] = {indices[3 * geomID + 0], 
                    indices[3 * geomID + 1],
                    indices[3 * geomID + 2]};

    v0[0] = vertices[3 * index[0] + 0];
    v0[1] = vertices[3 * index[0] + 1];
    v0[2] = vertices[3 * index[0] + 2];

    v1[0] = vertices[3 * index[1] + 0];
    v1[1] = vertices[3 * index[1] + 1];
    v1[2] = vertices[3 * index[1] + 2];

    v2[0] = vertices[3 * index[2] + 0];
    v2[1] = vertices[3 * index[2] + 1];
    v2[2] = vertices[3 * index[2] + 2];

    Real maxim[3], minimum[3];
    maxim[0] = std::max({v0[0], v1[0], v2[0]});
    maxim[1] = std::max({v0[1], v1[1], v2[1]});
    maxim[2] = std::max({v0[2], v1[2], v2[2]});
    minimum[0] = std::min({v0[0], v1[0], v2[0]});
    minimum[1] = std::min({v0[1], v1[1], v2[1]});
    minimum[2] = std::min({v0[2], v1[2], v2[2]});
    aabb.Setter(maxim, minimum);
}

template class SphereData<float>;
template class SphereData<double>;
template class TriangleData<float>;
template class TriangleData<double>;

template class BinaryBVH<float, SphereData<float>>;
template class BinaryBVH<double, SphereData<double>>;
template class BinaryBVH<float, TriangleData<float>>;
template class BinaryBVH<double, TriangleData<double>>;