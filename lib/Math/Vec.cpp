#include "Vec.hpp"

template<class Real>
Vec3<Real>::Vec3(Real x, Real y, Real z)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

template<class Real>
Vec3<Real>::Vec3(const Real* a)
{
    v[0] = a[0];
    v[1] = a[1];
    v[2] = a[2];
}


template<class Real>
Vec3<Real> Vec3<Real>::operator+(Vec3 vec)
{
    return Vec3(v[0] + vec[0],
                v[1] + vec[1],
                v[2] + vec[2]);
}

template<class Real>
Vec3<Real> Vec3<Real>::operator-(Vec3 vec)
{
    return Vec3(v[0] - vec[0],
                v[1] - vec[1],
                v[2] - vec[2]);
}

template<class Real>
Vec3<Real> Vec3<Real>::operator*(Vec3 vec)
{
    return Vec3(v[0] * vec[0],
                v[1] * vec[1],
                v[2] * vec[2]);
}

template<class Real>
Vec3<Real> operator*(Real t, Vec3<Real> vec)
{
    return Vec3(t * vec[0], 
                t * vec[1], 
                t * vec[2]);
}

template<class Real>
Vec3<Real> operator*(Vec3<Real> vec, Real t)
{
    return t * vec;
}

template<class Real>
Real Vec3<Real>::length()
{
    return std::sqrt(v[0] * v[0] + 
                     v[1] * v[1] +
                     v[2] * v[2]);
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
Vec3<Real> Vec3<Real>::normalize()
{
    return Vec3(v[0], v[1], v[2])
            * 1.0 / std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

template<class Real>
std::ostream &operator<<(std::ostream& stream, const Vec3<Real>& v)
{
    stream << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
    return stream;
}

template<class Real>
void ONB(const Vec3<Real>& n, Vec3<Real>& b1, Vec3<Real>& b2)
{
    if(n[1] < 0)
    {
        const Real a = 1.0/(1.0 - n[1]);
        const Real b = n[0] * n[2] * a;
        b1 = Vec3<Real>(1.0 - n[0] * n[0] * a, n[0], -b);
        b2 = Vec3<Real>(b, -n[2], n[2]*n[2]*a-1.0);
    }
    else
    {
        const Real a = 1.0/(1.0 + n[0]);
        const Real b = -n[0]*n[2]*a;
        b1 = Vec3<Real>(1.0 - n[0] * n[0] * a, -n[0], b);
        b2 = Vec3<Real>(b, -n[2], 1.0 - n[2]*n[2]*a);
    }
}

template<class Real>
Vec3<Real> worldtolocal(const Vec3<Real>& w, const Vec3<Real>& e0, const Vec3<Real>& e1, const Vec3<Real>& e2) {
    return Vec3(dot(e0, w), dot(e1, w), dot(e2, w));
}



template class Vec3<float>;
template class Vec3<double>;