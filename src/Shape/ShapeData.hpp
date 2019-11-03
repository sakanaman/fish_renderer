#ifndef SHAPEDATA_HPP
#define SHAPEDATA_HPP

#include <Accelerator.hpp>

template<class Real>
class TriangleData
{
public:
    bool intersect();
private:
    Real* vertices;
    Real* indices;
};

#endif