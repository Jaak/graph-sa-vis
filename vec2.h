#pragma once

#include <cmath>

#include <boost/optional.hpp>

// TODO: generalize

template <typename T>
struct vec2 {
    T x;
    T y;

    vec2(T x, T y)
        : x(x), y(y)
    { }
};

using vec2f = vec2<float>;

inline float sqr(float x) { return x*x; }
inline vec2f operator + (vec2f u, vec2f v) { return {u.x + v.x, u.y + v.y}; }
inline vec2f operator - (vec2f u) { return {- u.x, - u.y}; }
inline vec2f operator - (vec2f u, vec2f v) { return {u.x - v.x, u.y - v.y}; }
inline vec2f operator * (float s, vec2f v) { return {s * v.x, s * v.y}; }
inline vec2f operator * (vec2f u, float s) { return {u.x * s, u.y * s}; }
inline vec2f operator / (vec2f u, float s) { return {u.x / s, u.y / s}; }
inline float dot(vec2f u, vec2f v) { return u.x*v.x + u.y*v.y; }
inline float cross(vec2f u, vec2f v) { return u.x*v.y - u.y*v.x; }
inline float sqrlen(vec2f u) { return dot(u, u); }
inline float len(vec2f u) { return std::sqrt(sqrlen(u)); }
inline vec2f normalised(vec2f u) { return u / len(u); }
inline float sqrdist(vec2f u, vec2f v) { return sqrlen(u - v); }
inline float dist(vec2f u, vec2f v) { return std::sqrt(sqrdist(u, v)); }

// Checks if point q lies parallel to the line segment [p1, p2].
// If so, returns the square of its distance from the line.
inline boost::optional<float> sqrDistanceFromLine(vec2f p1, vec2f p2, vec2f q) {
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
    if (x == 0)
        return 0;
    return x > 0 ? 1 : 2;
}

#define ON_SEGMENT(p, q, r) (q.x <= std::fmax(p.x, r.x) && q.x >= std::fmin(p.x, r.x) && q.y <= std::fmax(p.y, r.y) && q.y >= std::fmin(p.y, r.y))
#define ORIENTATION(p, q, r) mysign((q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y))

inline bool intersects(vec2f p1, vec2f q1, vec2f p2, vec2f q2) {
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