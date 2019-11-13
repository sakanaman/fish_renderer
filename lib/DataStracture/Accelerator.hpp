#ifndef ACCEL_HPP
#define ACCEL_HPP

#include <memory>
#include <vector>


template<class Real>
class Hit
{
public:
    Hit();
    Real* GetPos();
    Real* GetNg();
    Real GetT() const;
    void SetT(const Real _t);
    int GetID() const;
    void SetPos(const Real* pos);
    void SetNg(const Real* n);
    void SetID(const int num);
private:
    Real point[3];
    Real normal[3];
    int geomID;
    Real t;
};

enum class Axis;
enum class LR;

template<class Real>
class AABB
{
public:
    AABB(const Real* max,const Real* min);
    AABB();
    void Setter(const Real* max, const Real* min);
    bool intersect(const Real* ray_origin, const Real* ray_dir) const;
    AABB Union(const AABB& r) const;
    AABB Union(const Real a[3]) const;
    Real GetMax(const Axis axis) const;
    Real GetMin(const Axis axis) const;
    Real AreaAABB() const;
private:
    // Real max[3] =  {std::numeric_limits<Real>::lowest(),
    //                 std::numeric_limits<Real>::lowest(),
    //                 std::numeric_limits<Real>::lowest()};;
    // Real min[3] =  {std::numeric_limits<Real>::max(),
    //                 std::numeric_limits<Real>::max(),
    //                 std::numeric_limits<Real>::max()};
    Real max[3] = {-1000000,
                   -1000000,
                   -1000000};
    Real min[3] = {1000000,
                   1000000,
                   1000000,};
};


template<class Real>
class BBVHNode
{
public:
    void SetSplitDimention(const Axis axis);
    void SetChildIndex(const LR lr,const int index);
    void SetLeafNode(const int* range);
    void SetAABB(const AABB<Real>& aabb);
    bool isLeaf() const;
    int GetBegin() const;
    int GetEnd() const;
    int GetChild(const LR lr);
    AABB<Real> GetAABB() const;
private:
    AABB<Real> aabb;
    int childinfo_left = 0;
    int childinfo_right = 0;
};

enum class Evaluator
{
    SAH
};
template<class Real, class ShapeData>
class BinaryBVH
{
public:
    BinaryBVH(const ShapeData& shapedata, const int face_num);
    int PartitionSAH(const int* range, const AABB<Real>& bigaabb);
    void BuildBVH(const Evaluator split_way);
    bool Traverse(Hit<Real>& hit, const Real* ray_origin, const Real* ray_dir, const int index) const;
private:
    ShapeData shapedata;
    std::unique_ptr<BBVHNode<Real>[]> bvh_nodes;
    std::unique_ptr<AABB<Real>[]> aabbs;
    std::vector<int> faces;
};

//shapedata

template<class Real>
class TriangleData
{
public:
    TriangleData();
    TriangleData(Real* vertices, int* indices);
    bool intersect(const int geomID, const Real* ray_origin, const Real* ray_dir, Hit<Real>& hit) const;
    void makeAABB(const int geomID, AABB<Real>& aabb) const;
private:
    Real* vertices;
    int* indices;
};

template<class Real>
class SphereData
{
public:
    SphereData();
    SphereData(Real* rad_cents);
    bool intersect(const int geomID, const Real* ray_origin , const Real* ray_dir, Hit<Real>& hit) const;
    void makeAABB(const int geomID, AABB<Real>& aabb) const;
private:
    Real* rad_cents;
    //__________________________________________
    //| radius | center.x | center.y | center.z | ...
    //------------------------------------------
    //                per geomID
};
#endif