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
Vec3<Real> Vec3<Real>::operator+(const Vec3& vec) const
{
    return Vec3(v[0] + vec[0],
                v[1] + vec[1],
                v[2] + vec[2]);
}

template<class Real>
Vec3<Real> Vec3<Real>::operator-(const Vec3& vec) const
{
    return Vec3(v[0] - vec[0],
                v[1] - vec[1],
                v[2] - vec[2]);
}

template<class Real>
Vec3<Real> Vec3<Real>::operator*(const Vec3& vec) const
{
    return Vec3(v[0] * vec[0],
                v[1] * vec[1],
                v[2] * vec[2]);
}

template<class Real>
Vec3<Real> Vec3<Real>::operator/(const Vec3& vec) const
{
    return Vec3(v[0] / vec[0],
                v[1] / vec[1],
                v[2] / vec[2]);
}

template<class Real>
Vec3<Real> operator*(const Real t, const Vec3<Real>& vec)
{
    return Vec3<Real>(t * vec[0], 
                      t * vec[1], 
                      t * vec[2]);
}

template<class Real>
Vec3<Real> operator*(const Vec3<Real>& vec, const Real t)
{
    return t * vec;
}

template<class Real>
Real Vec3<Real>::length() const
{
    return std::sqrt(v[0] * v[0] + 
                     v[1] * v[1] +
                     v[2] * v[2]);
}

template<class Real>
Real Vec3<Real>::dot(const Vec3& vec) const
{
    return v[0] * vec[0] + 
           v[1] * vec[1] + 
           v[2] * vec[2];
}

template<class Real>
Vec3<Real> Vec3<Real>::cross(const Vec3& vec) const
{
    return Vec3(v[1] * vec[2] - vec[1] * v[2],
                v[2] * vec[0] - vec[2] * v[0],
                v[0] * vec[1] - vec[0] * v[1]);
}

template<class Real>
Vec3<Real> Vec3<Real>::normalize() const
{
    return Vec3(v[0], v[1], v[2])
            * (static_cast<Real>(1.0) / std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
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
    Real sign = copysign(1.0, n[1]);
    const Real a = -1.0 / (sign + n[1]);
    const Real b = n[0] * n[2] * a;
    b1 = Vec3<Real>(1.0 + sign * n[0] * n[0] * a,  -sign * n[0] ,sign * b);
    b2 = Vec3<Real>(b, -n[2], sign + n[2] * n[2] * a);

    // Real sign = std::copysign(1.0, n[2]);
    // const Real a = -1.0 / (sign + n[2]);
    // const Real b = n[0] * n[1] * a;
    // b1 = Vec3<Real>(1.0 + sign * n[0] * n[0] * a, sign * b, -sign * n[0] );
    // b2 = Vec3<Real>(b, sign + n[1] * n[1] * a, -n[1]);

    // std::cout << "n's length: " << n.length() << std::endl;
    // std::cout << "b1's length: " << b1.length() << std::endl;
    // std::cout << "b2's length: " << b2.length() << std::endl;
}

template<class Real>
Vec3<Real> worldtolocal(const Vec3<Real>& w, const Vec3<Real>& e0, const Vec3<Real>& e1, const Vec3<Real>& e2) {
    return Vec3<Real>(e0.dot(w), e1.dot(w), e2.dot(w));
}


template class Vec3<float>;
template class Vec3<double>;
template Vec3<float> operator*(const float t, const Vec3<float>& vec);
template Vec3<float> operator*(const Vec3<float>& vec, const float t);
template std::ostream& operator<<(std::ostream&, const Vec3<float>&);
template std::ostream& operator<<(std::ostream&, const Vec3<double>&);
template void ONB(const Vec3<float>& n, Vec3<float>& b1, Vec3<float>& b2);
template void ONB(const Vec3<double>& n, Vec3<double>& b1, Vec3<double>& b2);
template Vec3<double> worldtolocal(const Vec3<double>& w, const Vec3<double>& e0, const Vec3<double>& e1, const Vec3<double>& e2);
template Vec3<float> worldtolocal(const Vec3<float>& w, const Vec3<float>& e0, const Vec3<float>& e1, const Vec3<float>& e2);