#include "Accelerator.hpp"

//Hit class implementation
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
    point = pos;
}

template<class Real>
void Hit<Real>::SetNg(Real* n)
{
    normal = n;
}

template<class Real>
void Hit<Real>::SetID(int num)
{
    geomID = num;
}


//Ray class implementation
template<class Real>
Ray<Real>::Ray(Real* _origin, Real* _dir)
{
    origin = _origin;
    dir = _dir;
}

template<class Real>
void Ray<Real>::SetDir(Real* _dir)
{
    dir = _dir;
}

template<class Real>
void Ray<Real>::SetOrigin(Real* _origin)
{
    origin = _origin;
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


//AABB class implementation
template<class Real>
AABB<Real>::AABB(Real* _max, Real* _min)
{
    max = _max;
    min = _min;
}

template<class Real>
bool AABB<Real>::intersect(const Ray<Real>& ray) const{
    float t_max = std::numeric_limits<Real>::max();
    float t_min = std::numeric_limits<Real>::lowest();

    for(int i = 0; i < 3; ++i) {
        float t1 = (min[i] - ray.origin[i])/ray.dir[i];
        float t2 = (max[i] - ray.origin[i])/ray.dir[i];
        
        float t_near = std::min(t1, t2);
        float t_far = std::max(t1, t2);
        
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
Real AABB<Real>::AreaAABB()
{
    Real result;
    result = 2*((max[0] - min[0])*(max[1] - min[1]))
            +2*((max[1] - min[1])*(max[2] - min[2]))
            +2*((max[2] - min[2])*(max[0] - min[0]));
    return result;
}

//BinaryBVH`s Node implementation
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
void BBVHNode<Real>::SetAABB(const AABB<Real>& _aabb)
{
    aabb = _aabb;
}


//BinaryBVH implementation
