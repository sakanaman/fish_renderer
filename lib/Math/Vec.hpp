#ifndef VEC_HPP
#define VEC_HPP

#include <cmath>
#include <iostream>

template<class Real>
class Vec3 {
public:
    Vec3(Real x, Real y, Real z);
    Vec3(Real v[3]);
    explicit operator Real*() const;
    Vec3 operator+(Vec3 vec);
    Vec3 operator-(Vec3 vec);
    Vec3 operator*(Vec3 vec);
    Real dot(Vec3 vec);
    Vec3 cross(Vec3 vec);
private:
    Real v[3];
};

template<class Real>
std::ostream &operator<<(std::ostream& stream, const Vec3<Real>& v);

#endif