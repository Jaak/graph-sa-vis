#pragma once

#include <cmath>
#include <boost/optional.hpp>
#include <type_traits>

template <
    typename T,
    typename = typename std::enable_if<std::is_floating_point<T>::value>::type
>
struct vec2 {
    T x;
    T y;

    vec2(T x, T y)
        : x(x), y(y)
    { }

    friend inline vec2 operator + (vec2 u, vec2 v) { return {u.x + v.x, u.y + v.y}; }
    friend inline vec2 operator - (vec2 u) { return {- u.x, - u.y}; }
    friend inline vec2 operator - (vec2 u, vec2 v) { return {u.x - v.x, u.y - v.y}; }
    friend inline vec2 operator * (T s, vec2 v) { return {s * v.x, s * v.y}; }
    friend inline vec2 operator * (vec2 u, T s) { return {u.x * s, u.y * s}; }
    friend inline vec2 operator / (vec2 u, T s) { return {u.x / s, u.y / s}; }
};

template <typename T>
inline float dot(vec2<T> u, vec2<T> v) { return u.x*v.x + u.y*v.y; }

template <typename T>
inline float cross(vec2<T> u, vec2<T> v) { return u.x*v.y - u.y*v.x; }

template <typename T>
inline float sqrlen(vec2<T> u) { return dot(u, u); }

template <typename T>
inline float len(vec2<T> u) { return std::sqrt(sqrlen(u)); }

template <typename T>
inline vec2<T> normalised(vec2<T> u) { return u / len(u); }

template <typename T>
inline float sqrdist(vec2<T> u, vec2<T> v) { return sqrlen(u - v); }

template <typename T>
inline float dist(vec2<T> u, vec2<T> v) { return std::sqrt(sqrdist(u, v)); }


// Checks if point q lies parallel to the line segment [p1, p2].
// If so, returns the square of its distance from the line.
template <typename T>
inline boost::optional<T> sqrDistanceFromLine(vec2<T> p1, vec2<T> p2, vec2<T> q) {
    const auto u = normalised(p2 - p1);
    
    auto v1 = q - p1;
    const auto lSqr = sqrlen(v1);
    if (lSqr == 0.0)
        return boost::none;

    v1 = v1 / std::sqrt(lSqr);
    const auto cosA = dot(u, v1);
    if (cosA < 0.0 || cosA > 1.0)
        return boost::none;

    const auto v2 = normalised(q - p2);
    const auto cosB = dot(-u, v2);
    if (cosB < 0.0 || cosB > 1.0)
        return boost::none;

    return lSqr * (1 - cosA*cosA);
}

inline int mysign(float x) {
    return x == 0 ? 0 : (x > 0 ? 1 : 2);
}

#define ON_SEGMENT(p, q, r) (q.x <= std::fmax(p.x, r.x) && q.x >= std::fmin(p.x, r.x) && q.y <= std::fmax(p.y, r.y) && q.y >= std::fmin(p.y, r.y))
#define ORIENTATION(p, q, r) mysign((q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y))

template <typename T>
inline bool intersects(vec2<T> p1, vec2<T> q1, vec2<T> p2, vec2<T> q2) {
    const auto o1 = ORIENTATION(p1, q1, p2);
    const auto o2 = ORIENTATION(p1, q1, q2);
    const auto o3 = ORIENTATION(p2, q2, p1);
    const auto o4 = ORIENTATION(p2, q2, q1);

    bool intersect = false;
    intersect |= (o1 != o2 && o3 != o4);
    intersect |= (o1 == 0 && ON_SEGMENT(p1, p2, q1));
    intersect |= (o2 == 0 && ON_SEGMENT(p1, q2, q1));
    intersect |= (o3 == 0 && ON_SEGMENT(p2, p1, q2));
    intersect |= (o4 == 0 && ON_SEGMENT(p2, q1, q2));
    return intersect;
}

#undef ON_SEGMENT
#undef ORIENTATION

using vec2f = vec2<float>;
using vec2d = vec2<double>;