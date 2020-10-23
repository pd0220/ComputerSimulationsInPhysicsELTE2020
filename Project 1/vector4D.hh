// 4D vector structure implementation for further use with
// second order differential equations

// used headers and libraries
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

// 4D vector struct
template <typename T>
struct vector4
{
    T x1, x2, v1, v2;

    // setting vector values (=)
    vector4<T> &operator=(vector4<T> const &v)
    {
        x1 = v.x1;
        x2 = v.x2;
        v1 = v.v1;
        v2 = v.v2;
        return *this;
    }

    // summation (+=)
    template <typename t>
    auto &operator+=(vector4<t> const &v)
    {
        x1 += (t)v.x1;
        x2 += (t)v.x2;
        v1 += (t)v.v1;
        v2 += (t)v.v2;
        return *this;
    }

    // substraction (-=)
    template <typename t>
    auto &operator-=(vector4<t> const &v)
    {
        x1 -= (t)v.x1;
        x2 -= (t)v.x2;
        v1 -= (t)v.v1;
        v2 -= (t)v.v2;
        return *this;
    }

    // multiplication by scalar (*=)
    template <typename t>
    auto &operator*=(t const &a)
    {
        x1 *= (t)a;
        x2 *= (t)a;
        v1 *= (t)a;
        v2 *= (t)a;
        return *this;
    }

    // division by scalar (/=)
    template <typename t>
    auto &operator/=(t const &a)
    {
        x1 /= (t)a;
        x2 /= (t)a;
        v1 /= (t)a;
        v2 /= (t)a;
        return *this;
    }
};

// output stream
template <typename T>
std::ostream &operator<<(std::ostream &o, vector4<T> const &v)
{
    o << v.x1 << ", " << v.x2 << ", " << v.v1 << ", " << v.v2;
    return o;
}

// summation of vectors (+)
template <typename T1, typename T2>
auto operator+(vector4<T1> const &v1, vector4<T2> const &v2)
{
    using R = decltype(v1.x1 + v2.x1);
    return vector4<R>{v1.x1 + v2.x1, v1.x2 + v2.x2, v1.v1 + v2.v1, v1.v2 + v2.v2};
}

// substraction of vectors (-)
template <typename T1, typename T2>
auto operator-(vector4<T1> const &v1, vector4<T2> const &v2)
{
    using R = decltype(v1.x1 - v2.x1);
    return vector4<R>{v1.x1 - v2.x1, v1.x2 - v2.x2, v1.v1 + v2.v1, v1.v2 + v2.v2};
}

// multiplication by a scalar (from right) (*)
template <typename Tv, typename T>
auto operator*(vector4<Tv> const &v, T const &a)
{
    using R = decltype(v.x1 * a);
    return vector4<R>{v.x1 * a, v.x2 * a, v.v1 * a, v.v2 * a};
}

// multiplication by a scalar (from left) (*)
template <typename Tv, typename T>
auto operator*(T const &a, vector4<Tv> const &v)
{
    using R = decltype(a * v.x1);
    return vector4<R>{a * v.x1, a * v.x2, a * v.v1, a * v.v2};
}

// division by a scalar (/)
template <typename Tv, typename T>
auto operator/(vector4<Tv> const &v, T const &a)
{
    using R = decltype(v.x1 / a);
    return vector4<R>{v.x1 / a, v.x2 / a, v.v1 / a, v.v2 / a};
}

// dot product
template <typename T1, typename T2>
auto dot(vector4<T1> const &v1, vector4<T2> const &v2)
{
    return v1.x1 * v2.x1 + v1.x2 * v1.x2 + v1.v1 * v2.v1 + v1.v2 * v2.v2;
}

// length of a vector
template <typename T>
auto length(vector4<T> const &v)
{
    return std::sqrt(v.x1 * v.x1 + v.x2 * v.x2 + v.v1 * v.v1 + v.v2 * v.v2);
}

// square length of a vector
template <typename T>
auto sqlength(vector4<T> const &v)
{
    return v.x1 * v.x1 + v.x2 * v.x2 + v.v1 * v.v1 + v.v2 * v.v2;
}

// normalize a vector
template <typename T>
vector4<T> normalize(vector4<T> const &v)
{
    if (v.x1 == 0 && v.x2 == 0 && v.v1 == 0 && v.v2 == 0)
    {
        std::cout << "error\nLenght of vector cannot be interpreted." << std::endl;
        exit(-1);
    }
    return v / length(v);
}