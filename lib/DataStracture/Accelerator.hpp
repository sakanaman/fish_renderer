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
    Real GetT();
    void SetT(Real _t);
    int GetID();
    void SetPos(Real* pos);
    void SetNg(Real* n);
    void SetID(int num);
private:
    Real point[3];
    Real normal[3];
    int geomID;
    Real t;
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

enum class Axis;
enum class LR;

template<class Real>
class AABB
{
public:
    AABB(Real* max, Real* min);
    AABB();
    void Setter(Real* max, Real* min);
    bool intersect(Ray<Real>& ray) const;
    AABB Union(const AABB& r);
    AABB Union(const Real a[3]);
    Real GetMax(Axis axis);
    Real GetMin(Axis axis);
    Real AreaAABB();
private:
    Real max[3] =  {std::numeric_limits<Real>::lowest(),
                    std::numeric_limits<Real>::lowest(),
                    std::numeric_limits<Real>::lowest()};;
    Real min[3] =  {std::numeric_limits<Real>::max(),
                    std::numeric_limits<Real>::max(),
                    std::numeric_limits<Real>::max()};
};


template<class Real>
class BBVHNode
{
public:
    void SetSplitDimention(Axis axis);
    void SetChildIndex(LR lr, int index);
    void SetLeafNode(int range[3]);
    void SetAABB(const AABB<Real>& aabb);
    bool isLeaf();
    int GetBegin();
    int GetEnd();
    int GetChild(LR lr);
    AABB<Real> GetAABB(); 
private:
    AABB<Real> aabb;
    int childinfo_left = 0;
    int childinfo_right = 0;
};

enum class Evaluator;
template<class Real, class ShapeData>
class BinaryBVH
{
public:
    BinaryBVH(const ShapeData& shapedata, int face_num);
    int PartitionSAH(int range[3], AABB<Real>& bigaabb);
    void BuildBVH(Evaluator split_way);
    bool Traverse(Hit<Real>& hit, Ray<Real>& ray, int index);
private:
    ShapeData shapedata;
    std::unique_ptr<BBVHNode<Real>[]> bvh_nodes;
    std::unique_ptr<AABB<Real>[]> aabbs;
    std::vector<int> faces;
};

//shapedata

// template<class Real>
// class TriangleData
// {
// public:
//     TriangleData(Real* vertices, Real* indices);
//     bool intersect(int geomID, const Ray<Real>& ray, Hit<Real>& hit);
//     void makeAABB(int geomID, AABB<Real>& aabb);
// private:
//     Real* vertices;
//     Real* indices;
// };

template<class Real>
class SphereData
{
public:
    SphereData();
    SphereData(Real* rad_cents);
    bool intersect(int geomID, Ray<Real>& ray, Hit<Real>& hit);
    void makeAABB(int geomID, AABB<Real>& aabb);
private:
    Real* rad_cents;
    //__________________________________________
    //| radius | center.x | center.y | center.z | ...
    //------------------------------------------
    //                per geomID
};
#endif