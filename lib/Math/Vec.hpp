#ifndef VEC_HPP
#define VEC_HPP

#include <cmath>
#include <iostream>

template<class Real>
class Vec3 {
public:
    Vec3(Real x, Real y, Real z);
    Vec3(const Real* v);
    Vec3(){};
    Vec3 operator+(const Vec3& vec) const;
    Vec3 operator-(const Vec3& vec) const;
    Vec3 operator*(const Vec3& vec) const;
    Vec3 operator/(const Vec3& vec) const;
    Real length() const;
    Real dot(const Vec3& vec) const;
    Vec3 cross(const Vec3& vec) const;
    Vec3 normalize() const;
    Real& operator[](int i) {return v[i];};
    Real operator[](int i) const {return v[i];};
private:
    Real v[3];
};

template<class Real>
Vec3<Real> operator*(const Real t, const Vec3<Real>& vec);

template<class Real>
Vec3<Real> operator*(const Vec3<Real>& vec, const Real t);

template<class Real>
std::ostream &operator<<(std::ostream& stream, const Vec3<Real>& v);

template<class Real>
void ONB(const Vec3<Real>& n, Vec3<Real>& b1, Vec3<Real>& b2);

template<class Real>
Vec3<Real> worldtolocal(const Vec3<Real>& w, const Vec3<Real>& e0, const Vec3<Real>& e1, const Vec3<Real>& e2);

#endif