#include "Accelerator.hpp"
#include <cmath>

//Hit class implementation
template<class Real>
Real Hit<Real>::GetT()
{
    return t;
}

template<class Real>
void Hit<Real>::SetT(Real _t)
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
int Hit<Real>::GetID()
{
    return geomID;
}

template<class Real>
void Hit<Real>::SetPos(Real* pos)
{   
    point[0] = pos[0];
    point[1] = pos[1];
    point[2] = pos[2];
}

template<class Real>
void Hit<Real>::SetNg(Real* n)
{
    normal[0] = n[0];
    normal[1] = n[1];
    normal[2] = n[2];
}

template<class Real>
void Hit<Real>::SetID(int num)
{
    geomID = num;
}

template class Hit<float>;
template class Hit<double>;

//Ray class implementation
template<class Real>
Ray<Real>::Ray(Real* _origin, Real* _dir)
{
    origin[0] = _origin[0];
    origin[1] = _origin[1];
    origin[2] = _origin[2];

    dir[0] = _dir[0];
    dir[1] = _dir[1];
    dir[2] = _dir[2];
}

template<class Real>
void Ray<Real>::SetDir(Real* _dir)
{
    dir[0] = _dir[0];
    dir[1] = _dir[1];
    dir[2] = _dir[2];
}

template<class Real>
void Ray<Real>::SetOrigin(Real* _origin)
{
    origin[0] = _origin[0];
    origin[1] = _origin[1];
    origin[2] = _origin[2];
}

template<class Real>
Real* Ray<Real>::GetDir()
{
    return dir;
}

template<class Real>
Real* Ray<Real>::GetOrigin()
{
    return origin;
}

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
AABB<Real>::AABB(Real* _max, Real* _min)
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
void AABB<Real>::Setter(Real* _max, Real* _min)
{
    max[0] = _max[0];
    max[1] = _max[1];
    max[2] = _max[2];

    min[0] = _min[0];
    min[1] = _min[1];
    min[2] = _min[2];
}

template<class Real>
bool AABB<Real>::intersect(Ray<Real>& ray) const {
    Real t_max = std::numeric_limits<Real>::max();
    Real t_min = std::numeric_limits<Real>::lowest();

    for(int i = 0; i < 3; ++i) {
        Real t1 = (min[i] - (ray.GetOrigin())[i])/(ray.GetDir())[i];
        Real t2 = (max[i] - (ray.GetOrigin())[i])/(ray.GetDir())[i];
        
        Real t_near = std::min(t1, t2);
        Real t_far = std::max(t1, t2);
        
        t_max = std::min(t_max, t_far);
        t_min = std::max(t_min, t_near);
    
        if (t_min > t_max) return false;
    }
    return true;
}

template<class Real>
AABB<Real> AABB<Real>::Union(const AABB& r)
{
    AABB result;

    for(int i = 0; i < 3; ++i) {
        result.max[i] =  std::max(max[i], r.max[i]);
        result.min[i] =  std::min(min[i], r.min[i]);
    }

    return result;
}

template<class Real>
AABB<Real> AABB<Real>::Union(const Real a[3])
{
    AABB result;

    for(int i = 0; i < 3; ++i) {
        result.max[i] = std::max(max[i], a[i]);
        result.min[i] = std::min(min[i], a[i]);
    }
    return result;
}

template<class Real>
Real AABB<Real>::GetMax(Axis axis)
{
    return max[static_cast<int>(axis)];
}

template<class Real>
Real AABB<Real>::GetMin(Axis axis)
{
    return min[static_cast<int>(axis)];
}

template<class Real>
Real AABB<Real>::AreaAABB()
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
void BBVHNode<Real>::SetSplitDimention(Axis axis)
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
void BBVHNode<Real>::SetChildIndex(LR lr, int index)
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
void BBVHNode<Real>::SetLeafNode(int range[3])
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
AABB<Real> BBVHNode<Real>::GetAABB()
{
    return aabb;
}

template<class Real>
int BBVHNode<Real>::GetBegin()
{
    return childinfo_left;
}

template<class Real>
int BBVHNode<Real>::GetEnd()
{
    return childinfo_right;
}

template<class Real>
bool BBVHNode<Real>::isLeaf()
{
    return !(childinfo_left & (1 << 31));
}

template<class Real>
int BBVHNode<Real>::GetChild(LR lr)
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
enum class Evaluator
{
    SAH
};

template<class Real, class ShapeData>
BinaryBVH<Real, ShapeData>::BinaryBVH(const ShapeData& _shapedata, int face_num)
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
int BinaryBVH<Real, ShapeData>::PartitionSAH(int range[3], AABB<Real>& bigaabb)
{
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
                centroid[j] = aabbs[faces[i]].GetMax(static_cast<Axis>(j));
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
                    // buckets[b].aabb.SetAABB((buckets[b].aabb.GetAABB()).Union(aabbs[faces[i]]));
                    buckets[b].aabb = (buckets[b].aabb).Union(aabbs[faces[i]]);
                }

                Real cost[nBuckets -1];
                for(int i = 0; i < nBuckets -1; ++i)
                {
                    AABB<Real> b0, b1;
                    int count0 = 0, count1 = 0;
                    for(int j = 0; j <= i; j++)
                    {
                        b0.Union(buckets[j].aabb);
                        count0 += buckets[j].count;
                    }
                    for(int j = i+1; j < nBuckets; j++)
                    {
                        b1.Union(buckets[j].aabb);
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
                                                        ((aabbs[faces[geomID]].GetMax(static_cast<Axis>(dim)) + aabbs[faces[geomID]].GetMin(static_cast<Axis>(dim))) * 0.5 + centoroidAABB.GetMin(static_cast<Axis>(dim)))
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
void BinaryBVH<Real, ShapeData>::BuildBVH(Evaluator eval)
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
                return;
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
        }
    }
}

template<class Real, class ShapeData>
bool BinaryBVH<Real, ShapeData>::Traverse(Hit<Real>& hit, Ray<Real>& ray, int index)
{
    bool is_hit_aabb = (bvh_nodes[index].GetAABB()).intersect(ray);

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
                if(shapedata.intersect(faces[i], ray, hit_each))
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
            bool is_hit1 = Traverse(hit, ray, bvh_nodes[index].GetChild(LR::Left));
            bool is_hit2 = Traverse(hit, ray, bvh_nodes[index].GetChild(LR::Right));
            return is_hit1 || is_hit2;
        }
        
    }
    
}


//shapedata implementation
template<class Real>
SphereData<Real>::SphereData(Real* rad_cents):rad_cents(rad_cents){}

template<class Real>
SphereData<Real>::SphereData(){}

template<class Real>
bool SphereData<Real>::intersect(int geomID, Ray<Real>& ray, Hit<Real>& hit) //reference from smallpt
{
    Real radius   = rad_cents[4 * geomID + 0];
    Real center_x = rad_cents[4 * geomID + 1];
    Real center_y = rad_cents[4 * geomID + 2];
    Real center_z = rad_cents[4 * geomID + 3];

    Real* ray_dir = ray.GetDir();
    Real* ray_origin = ray.GetOrigin();

    Real op[3] = {center_x - ray_origin[0],
                  center_y - ray_origin[1],
                  center_z - ray_origin[2]};

    Real t, eps = 1e-4;
    Real b = op[0] * ray_dir[0]
            +op[1] * ray_dir[1]
            +op[2] * ray_dir[2];
    
    Real det = b * b - (  op[0] * op[0] 
                        + op[1] * op[1]
                        + op[2] * op[2]) + radius * radius;
    if(det < 0) return false;
    else det = std::sqrt(det);
    Real hit_distance = (t = b - det)>eps ? t : ((t = b + det)>eps ? t : 0.0);
    if(hit_distance == 0.0) return false;
    Real hitpos[3] = {hit_distance * ray_dir[0] + ray_origin[0],
                      hit_distance * ray_dir[1] + ray_origin[1],
                      hit_distance * ray_dir[2] + ray_origin[2]};
    hit.SetPos(hitpos);
    Real hitnormal[3] = {(hitpos[0] - center_x)/radius, (hitpos[1] - center_y)/radius, (hitpos[2] - center_z)/radius}; 
    hit.SetNg(hitnormal);
    hit.SetID(geomID);
    hit.SetT(hit_distance);
    return true;
}

template<class Real>
void SphereData<Real>::makeAABB(int geomID, AABB<Real>& aabb)
{
    Real radius   = rad_cents[4 * geomID + 0];
    Real center_x = rad_cents[4 * geomID + 1];
    Real center_y = rad_cents[4 * geomID + 2];
    Real center_z = rad_cents[4 * geomID + 3];

    Real max[] = {center_x + radius, center_y + radius, center_z + radius};
    Real min[] = {center_x - radius, center_y - radius, center_z - radius};
    aabb.Setter(max, min);
}

template class SphereData<float>;
template class SphereData<double>;
template class BinaryBVH<float, SphereData<float>>;
template class BinaryBVH<double, SphereData<double>>;
