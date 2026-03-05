#pragma once
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>

// =============================================================================
//  geometry.h  —  Personal C++ Geometry Template
//  Author  : MI
//  Created : 2026-03
//  Repo    : https://github.com/LLLY248/geometry-template
//
//  更新记录：
//  2026-03  Part 0-3: 常量、Vec2、Vec3、基本几何判断
// =============================================================================


// =============================================================================
//  Part 0 : 常量与工具函数
// =============================================================================

namespace Geo {

constexpr float EPS = 1e-6f;	 // 浮点判零阈值，几何判断用
constexpr float PI = 3.14159265358979323846f;

// 浮点判零
inline bool isZero(float x) { return std::abs(x) < EPS; }

// 浮点判等
inline bool equal(float a, float b) { return std::abs(a - b) < EPS; }

// 符号函数：返回 -1 / 0 / 1
inline int sign(float x) {
    if (x > EPS) return  1;
    if (x < -EPS) return -1;
    return 0;
}


// =============================================================================
//  Part 1 : Vec2 — 二维向量
// =============================================================================

struct Vec2 {
    float x, y;

    // ---------- 构造 ----------
    Vec2() : x(0), y(0) {}
    Vec2(float x, float y) : x(x), y(y) {}

    // ---------- 基本运算 ----------
    Vec2 operator+(const Vec2& v) const { return { x + v.x, y + v.y }; }
    Vec2 operator-(const Vec2& v) const { return { x - v.x, y - v.y }; }
    Vec2 operator*(float f)       const { return { x * f,   y * f }; }
    Vec2 operator/(float f)       const { return { x / f,   y / f }; }

    Vec2& operator+=(const Vec2& v) { x += v.x; y += v.y; return *this; }
    Vec2& operator-=(const Vec2& v) { x -= v.x; y -= v.y; return *this; }
    Vec2& operator*=(float f) { x *= f;   y *= f;   return *this; }

    bool operator==(const Vec2& v) const { return equal(x, v.x) && equal(y, v.y); }
    bool operator!=(const Vec2& v) const { return !(*this == v); }

    // ---------- 向量运算 ----------

    // 点积
    float dot(const Vec2& v) const { return x * v.x + y * v.y; }

    // 2D叉积（标量）
    // 正值 → v 在 *this 的逆时针方向
    // 负值 → v 在 *this 的顺时针方向
    // 零   → 平行或共线
    float cross(const Vec2& v) const { return x * v.y - y * v.x; }

    // 向量长度
    float length()  const { return std::sqrt(x * x + y * y); }
    float length2() const { return x * x + y * y; }   // 长度的平方，避免开方更快

    // 单位向量（返回新向量，不修改自身）
    // 零向量时返回自身，不崩溃
    Vec2 normalized() const {
        float len = length();
        if (len < EPS) return *this;
        return *this / len;
    }

    // 逆时针旋转 90 度的垂直向量
    Vec2 perp() const { return { -y, x }; }

    // ---------- 工具 ----------
    bool isZero() const { return length2() < EPS * EPS; }
};

// 支持 scalar * Vec2 写法
inline Vec2 operator*(float f, const Vec2& v) { return v * f; }

// 两点距离
inline float dist(const Vec2& a, const Vec2& b) { return (a - b).length(); }
inline float dist2(const Vec2& a, const Vec2& b) { return (a - b).length2(); }


// =============================================================================
//  Part 2 : Vec3 — 三维向量
// =============================================================================

struct Vec3 {
    float x, y, z;

    // ---------- 构造 ----------
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x, float y, float z) : x(x), y(y), z(z) {}

    // ---------- 基本运算 ----------
    Vec3 operator+(const Vec3& v) const { return { x + v.x, y + v.y, z + v.z }; }
    Vec3 operator-(const Vec3& v) const { return { x - v.x, y - v.y, z - v.z }; }
    Vec3 operator*(float f)       const { return { x * f,   y * f,   z * f }; }
    Vec3 operator/(float f)       const { return { x / f,   y / f,   z / f }; }

    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    Vec3& operator*=(float f) { x *= f;   y *= f;   z *= f;   return *this; }

    bool operator==(const Vec3& v) const { return equal(x, v.x) && equal(y, v.y) && equal(z, v.z); }
    bool operator!=(const Vec3& v) const { return !(*this == v); }

    // ---------- 向量运算 ----------

    // 点积
    float dot(const Vec3& v) const {
        return x * v.x + y * v.y + z * v.z;
    }

    // 叉积（返回垂直于两向量的新向量）
    // 记忆方法：
    //   i 分量：y*vz - z*vy
    //   j 分量：z*vx - x*vz  （注意符号）
    //   k 分量：x*vy - y*vx
    Vec3 cross(const Vec3& v) const {
        return {
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        };
    }

    // 向量长度
    float length()  const { return std::sqrt(x * x + y * y + z * z); }
    float length2() const { return x * x + y * y + z * z; }

    // 单位向量（零向量时返回自身）
    Vec3 normalized() const {
        float len = length();
        if (len < EPS) return *this;
        return *this / len;
    }

    // ---------- 工具 ----------
    bool isZero() const { return length2() < EPS * EPS; }

    // 投影到另一个向量上的标量分量
    float projectOnto(const Vec3& axis) const {
        return this->dot(axis.normalized());
    }
};

inline Vec3 operator*(float f, const Vec3& v) { return v * f; }

inline float dist(const Vec3& a, const Vec3& b) { return (a - b).length(); }
inline float dist2(const Vec3& a, const Vec3& b) { return (a - b).length2(); }

// =============================================================================
//  Part 3 : 基本几何判断（2D）
//
//  所有函数都基于叉积，不依赖 atan2 或三角函数，精度更高速度更快
// =============================================================================

// -----------------------------------------------------------------------------
//  3.1  方向判断
// -----------------------------------------------------------------------------

// cross2D：以 O 为原点，计算向量 OA × OB
// 返回正值 → B 在 OA 的左侧（逆时针）
// 返回负值 → B 在 OA 的右侧（顺时针）
// 返回零   → O、A、B 三点共线
inline float cross2D(const Vec2& O, const Vec2& A, const Vec2& B) {
    return (A - O).cross(B - O);
}

// 判断 P 是否严格在有向线段 AB 的左侧
inline bool onLeft(const Vec2& A, const Vec2& B, const Vec2& P) {
    return cross2D(A, B, P) > EPS;
}

// 判断 P 是否严格在有向线段 AB 的右侧
inline bool onRight(const Vec2& A, const Vec2& B, const Vec2& P) {
    return cross2D(A, B, P) < -EPS;
}

// 判断三点是否共线
inline bool collinear(const Vec2& A, const Vec2& B, const Vec2& C) {
    return isZero(cross2D(A, B, C));
}

// -----------------------------------------------------------------------------
//  3.2  线段相关
// -----------------------------------------------------------------------------

// 判断点 P 是否在线段 AB 上（包括端点）
// 先用叉积判共线，再用点积判是否在端点之间
inline bool onSegment(const Vec2& A, const Vec2& B, const Vec2& P) {
    if (!isZero(cross2D(A, B, P))) return false;          // 不共线
    return (P - A).dot(P - B) <= EPS;                     // 在端点之间（含端点）
}

// 判断两线段 AB 和 CD 是否相交（包括端点接触）
inline bool segmentsIntersect(const Vec2& A, const Vec2& B,
    const Vec2& C, const Vec2& D) {
    float d1 = cross2D(C, D, A);
    float d2 = cross2D(C, D, B);
    float d3 = cross2D(A, B, C);
    float d4 = cross2D(A, B, D);

    // 一般情况：两段互相跨越
    if (sign(d1) * sign(d2) < 0 && sign(d3) * sign(d4) < 0) return true;

    // 退化情况：端点恰好落在另一条线段上
    if (onSegment(C, D, A)) return true;
    if (onSegment(C, D, B)) return true;
    if (onSegment(A, B, C)) return true;
    if (onSegment(A, B, D)) return true;

    return false;
}

// -----------------------------------------------------------------------------
//  3.3  三角形（2D）
// -----------------------------------------------------------------------------

// 三角形有符号面积（顶点逆时针时为正）
// 公式：叉积 / 2
inline float triangleSignedArea(const Vec2& A, const Vec2& B, const Vec2& C) {
    return cross2D(A, B, C) * 0.5f;
}

// 三角形面积（始终为正）
inline float triangleArea(const Vec2& A, const Vec2& B, const Vec2& C) {
    return std::abs(triangleSignedArea(A, B, C));
}

// 判断点 P 是否在三角形 ABC 内（包括边上）
// 原理：P 在 ABC 内 ⟺ P 在 AB、BC、CA 三条有向边的同侧
// 注意：无论 ABC 顶点是顺时针还是逆时针，has_neg && has_pos 的判断都正确
inline bool inTriangle(const Vec2& A, const Vec2& B, const Vec2& C, const Vec2& P) {
    float d1 = cross2D(A, B, P);
    float d2 = cross2D(B, C, P);
    float d3 = cross2D(C, A, P);
    bool has_neg = (d1 < -EPS) || (d2 < -EPS) || (d3 < -EPS);
    bool has_pos = (d1 > EPS) || (d2 > EPS) || (d3 > EPS);
    return !(has_neg && has_pos);
}

// 三角形外接圆圆心
// 原理：外接圆圆心是三边垂直平分线的交点
// 用途：Delaunay 三角剖分的核心计算
inline Vec2 circumcenter(const Vec2& A, const Vec2& B, const Vec2& C) {
    float ax = B.x - A.x, ay = B.y - A.y;
    float bx = C.x - A.x, by = C.y - A.y;
    float D = 2.0f * (ax * by - ay * bx);
    // D 接近零说明三点共线，外接圆不存在
    assert(!isZero(D) && "circumcenter: three points are collinear");
    float ux = (by * (ax * ax + ay * ay) - ay * (bx * bx + by * by)) / D;
    float uy = (ax * (bx * bx + by * by) - bx * (ax * ax + ay * ay)) / D;
    return { A.x + ux, A.y + uy };
}

// 判断点 P 是否在三角形 ABC 的外接圆内（严格内部）
// 用途：Delaunay 合法性检验（inCircle 判断）
// 原理：计算 3×3 行列式，正值表示 P 在外接圆内
// 前提：A、B、C 按逆时针顺序排列（否则结果取反）
inline bool inCircumcircle(const Vec2& A, const Vec2& B, const Vec2& C, const Vec2& P) {
    float ax = A.x - P.x, ay = A.y - P.y;
    float bx = B.x - P.x, by = B.y - P.y;
    float cx = C.x - P.x, cy = C.y - P.y;
    float det = ax * (by * (cx * cx + cy * cy) - cy * (bx * bx + by * by))
        - ay * (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by))
        + (ax * ax + ay * ay) * (bx * cy - by * cx);
    return det > EPS;
}

// -----------------------------------------------------------------------------
//  3.4  多边形（2D）
// -----------------------------------------------------------------------------

// 多边形有符号面积（顶点逆时针时为正）
// 公式：shoelace（鞋带公式）
inline float polygonSignedArea(const std::vector<Vec2>& poly) {
    float area = 0;
    int n = (int)poly.size();
    for (int i = 0; i < n; i++) {
        const Vec2& cur = poly[i];
        const Vec2& next = poly[(i + 1) % n];
        area += cur.cross(next);
    }
    return area * 0.5f;
}

inline float polygonArea(const std::vector<Vec2>& poly) {
    return std::abs(polygonSignedArea(poly));
}

// 射线法：判断点 P 是否在多边形内（含边上）
// 原理：向右发射水平射线，统计与多边形边的交叉次数
// 奇数次 → 在内部；偶数次 → 在外部
inline bool inPolygon(const std::vector<Vec2>& poly, const Vec2& P) {
    int n = (int)poly.size();
    if (n < 3) return false;

    // 先判断是否在某条边上
    for (int i = 0, j = n - 1; i < n; j = i++)
        if (onSegment(poly[i], poly[j], P)) return true;

    // 射线法
    bool inside = false;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        float yi = poly[i].y, yj = poly[j].y;
        float xi = poly[i].x, xj = poly[j].x;
        if ((yi > P.y) != (yj > P.y)) {
            float xIntersect = (xj - xi) * (P.y - yi) / (yj - yi) + xi;
            if (P.x < xIntersect)
                inside = !inside;
        }
    }
    return inside;
}

// -----------------------------------------------------------------------------
//  3.5  凸包（Graham Scan）
// -----------------------------------------------------------------------------

// 输入：任意一组 2D 点
// 输出：凸包顶点，按逆时针顺序排列
// 共线点处理：不保留共线中间点（<= 0 时弹出，等号处理了共线）
inline std::vector<Vec2> grahamScan(const std::vector<Vec2>& pts) {
    int n = (int)pts.size();
    if (n < 3) return pts;

    // Step 1：拷贝，找最低点（y 最小，相同则 x 最小）作为极角排序基准
    std::vector<Vec2> points(pts);
    int pivot = 0;
    for (int i = 1; i < n; i++)
        if (points[i].y < points[pivot].y ||
            (equal(points[i].y, points[pivot].y) && points[i].x < points[pivot].x))
            pivot = i;
    std::swap(points[0], points[pivot]);
    Vec2 p0 = points[0];

    // Step 2：按极角排序
    // 叉积 > 0 → a 极角更小 → a 排在前面
    // 共线时距离近的排前面
    std::sort(points.begin() + 1, points.end(), [&](const Vec2& a, const Vec2& b) {
        float c = cross2D(p0, a, b);
        if (!isZero(c)) return c > 0;
        return dist2(p0, a) < dist2(p0, b);
        });

    // Step 3：栈扫描
    // 叉积 <= 0 说明出现右转或共线，弹出栈顶
    std::vector<Vec2> hull;
    for (const auto& p : points) {
        while (hull.size() >= 2) {
            Vec2 top = hull.back();
            Vec2 second = hull[hull.size() - 2];
            if (cross2D(second, top, p) <= 0)
                hull.pop_back();
            else
                break;
        }
        hull.push_back(p);
    }
    return hull;
}

}   // namespace Geo