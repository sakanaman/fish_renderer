#ifndef ACCEL_HPP
#define ACCEL_HPP

#include <memory>
#include <vector>


template<class Real>
class Hit
{
public:
    Real* GetPos();
    Real* GetNg();
    int GetID();
    void SetPos(Real* pos);
    void SetNg(Real* n);
    void SetID(int num);
private:
    Real point[3];
    Real normal[3];
    int geomID;
};

template<class Real>
class Ray
{
public:
    Ray(Real* _origin, Real* _dir);
    void SetDir(Real* _dir);
    void SetOrigin(Real* _origin);
    Real* GetDir();
    Real* GetOrigin();
private:
    Real origin[3];
    Real dir[3];
};

template<class Real>
class AABB
{
public:
    AABB(Real* max, Real* min);
    bool intersect(const Ray<Real>& ray) const;
    AABB Union(const AABB& r);
    AABB Union(const Real a[3]);
    Real AreaAABB();
private:
    Real max[3] =  {std::numeric_limits<Real>::lowest(),
                    std::numeric_limits<Real>::lowest(),
                    std::numeric_limits<Real>::lowest()};;
    Real min[3] =  {std::numeric_limits<Real>::max(),
                    std::numeric_limits<Real>::max(),
                    std::numeric_limits<Real>::max()};
};

enum class Axis;
enum class LR;
template<class Real>
class BBVHNode
{
public:
    void SetSplitDimention(Axis axis);
    void SetChildIndex(LR lr, int index);
    void SetAABB(const AABB<Real>& aabb);
private:
    AABB<Real> aabb;
    int childinfo_left = 0;
    int childinfo_right = 0:
};

enum class Evaluator;
template<class Real, class ShapeData>
class BinaryBVH
{
public:
    BinaryBVH(ShapeData shapedata, int face_num);
    void BuildBVH(Evaluator split_way);
    bool Traverse(Hit<Real>& hit, const Ray<Real>& ray) const;
private:
    ShapeData shapedata;
    std::unique_ptr<BBVHNode<Real>[]> bvh_nodes;
    std::unique_ptr<AABB<Real>[]> aabbs;
    std::vector<int> faces;
};

#endif