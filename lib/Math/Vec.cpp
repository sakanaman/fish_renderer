#include "Vec.hpp"

template<class Real>
Vec3<Real>::Vec3(Real x, Real y, Real z)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

template<class Real>
Vec3<Real>::Vec3(Real a[3])
{
    v[0] = a[0];
    v[1] = a[1];
    v[2] = a[2];
}

template<class Real>
Vec3<Real>::operator Real *() const{
    return v;
}

template<class Real>
Vec3<Real> Vec3<Real>::operator+(Vec3 vec)
{
    return Vec3(v[0] + vec[0],
                v[1] + vec[1],
                v[2] + vec[2],);
}

template<class Real>
Vec3<Real> Vec3<Real>::operator-(Vec3 vec)
{
    return Vec3(v[0] - vec[0],
                v[1] - vec[1],
                v[2] - vec[2],);
}

template<class Real>
Vec3<Real> Vec3<Real>::operator*(Vec3 vec)
{
    return Vec3(v[0] * vec[0],
                v[1] * vec[1],
                v[2] * vec[2],);
}

template<class Real>
Real Vec3<Real>::dot(Vec3 vec)
{
    return v[0] * vec[0] + 
           v[1] * vec[1] + 
           v[2] * vec[2];
}

template<class Real>
Vec3<Real> Vec3<Real>::cross(Vec3 vec)
{
    return Vec3(v[1] * vec[2] - vec[1] * v[2],
                v[2] * vec[0] - vec[2] * v[0],
                v[0] * vec[1] - vec[0] * v[1]);
}

template<class Real>
std::ostream &operator<<(std::ostream& stream, const Vec3<Real>& v)
{
    stream << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
    return stream;
}

template class Vec3<float>;
template class Vec3<double>;