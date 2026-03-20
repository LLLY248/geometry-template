#pragma once
#include <cmath>
#include <cassert>
#include <vector>
#include <array>
#include <algorithm>
#include <map>
#include <fstream>
#include <string>
#include <limits>
// =============================================================================
//  geometry.h  —  Personal C++ Geometry Template
//  Author  : MI
//  Created : 2026-03
//  Repo    : https://github.com/LLLY248/geometry-template
//
//  更新记录：
//  2026-03  Part 0-3: 常量、Vec2、Vec3、基本几何判断
//  2026-03  Part 4:   Triangle3D
//  2026-03  Part 5:   HalfEdge 半边数据结构
//  2026-03  Part 6:   点定位
//  2026-03  Part 7:   Delaunay 三角剖分（Bowyer-Watson 完整实现）
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

    // 判断点 P 是否在三角形 ABC 的外接圆内（严格内部返回 true，圆上返回 false）
    //
    // 原理：将 A、B、C 平移使 P 成为原点，构造如下 3×3 行列式：
    //   | ax  ay  ax²+ay² |
    //   | bx  by  bx²+by² |
    //   | cx  cy  cx²+cy² |
    // 等价于把三点提升到抛物面 z=x²+y² 后判断四面体方向：
    //   det > 0  →  A B C 逆时针排列时，P 在外接圆内
    //   det < 0  →  P 在外接圆外
    //   det ≈ 0  →  四点共圆（退化情况）
    //
    // 前提：A B C 必须按逆时针排列，否则 det 符号取反
    // 用途：Delaunay 合法性检验（Lawson flip 的判断依据）
    inline bool inCircumcircle(const Vec2& A, const Vec2& B, const Vec2& C, const Vec2& P) {
        // 退化情况：三点共线，外接圆不存在，直接返回 false
        float area2 = cross2D(A, B, C);
        if (isZero(area2)) return false;

        // 如果顶点是顺时针排列，交换 B C 使其变为逆时针
        // 保证 det 的符号含义一致
        const Vec2* a = &A;
        const Vec2* b = &B;
        const Vec2* c = &C;
        if (area2 < 0) std::swap(b, c);

        // 平移：以 P 为原点，A B C 坐标相应减去 P
        float ax = a->x - P.x, ay = a->y - P.y;
        float bx = b->x - P.x, by = b->y - P.y;
        float cx = c->x - P.x, cy = c->y - P.y;

        // 各点到 P 的距离平方（提升到抛物面的 z 坐标）
        float ar2 = ax * ax + ay * ay;
        float br2 = bx * bx + by * by;
        float cr2 = cx * cx + cy * cy;

        // 3×3 行列式展开（按第三列余子式）
        float det = ax * (by * cr2 - cy * br2)
            - ay * (bx * cr2 - cx * br2)
            + ar2 * (bx * cy - by * cx);

        // 四点共圆时视为在圆外（不触发 flip），避免数值振荡
        if (isZero(det)) return false;

        return det > 0;       // 如果 det > 0，说明 P' 落在平面下方（对应圆内）
    }


    // -----------------------------------------------------------------------------
    //  3.3b Triangle2D — 二维三角形（对象封装版）
    //
    //  和 Part 4 的 Triangle3D 对称，提供一致的调用方式
    //  散点函数（triangleArea/inTriangle 等）仍然保留，Delaunay 算法里更常用
    // -----------------------------------------------------------------------------

    struct Triangle2D {
        Vec2 a, b, c;   // 三个顶点，建议逆时针排列

        // ---------- 构造 ----------
        Triangle2D() {}
        Triangle2D(const Vec2& a, const Vec2& b, const Vec2& c) : a(a), b(b), c(c) {}

        // ---------- 基本属性 ----------

        // 边向量（从顶点 a 出发）
        Vec2 edge_ab() const { return b - a; }
        Vec2 edge_ac() const { return c - a; }

        // 有符号面积（逆时针为正）
        float signedArea() const { return triangleSignedArea(a, b, c); }

        // 面积（始终为正）
        float area() const { return triangleArea(a, b, c); }

        // 重心
        Vec2 centroid() const { return (a + b + c) * (1.0f / 3.0f); }

        // 外接圆圆心
        Vec2 circumcenter() const { return Geo::circumcenter(a, b, c); }

        // ---------- 判断 ----------

        // 判断点 P 是否在三角形内（含边界）
        bool contains(const Vec2& P) const { return Geo::inTriangle(a, b, c, P); }

        // 判断点 P 是否在外接圆内
        bool inCircumcircle(const Vec2& P) const { return Geo::inCircumcircle(a, b, c, P); }

        // 顶点是否逆时针排列
        bool isCCW() const { return signedArea() > EPS; }
    };


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


    // =============================================================================
    //  Part 4 : Triangle3D — 三维三角形
    //
    //  重要原则：所有几何运算用边向量 (B-A), (C-A)，不用裸顶点坐标
    // =============================================================================

    struct Triangle3D {
        Vec3 a, b, c;   // 三个顶点，建议逆时针朝向正面

        // ---------- 构造 ----------
        Triangle3D() {}
        Triangle3D(const Vec3& a, const Vec3& b, const Vec3& c) : a(a), b(b), c(c) {}

        // ---------- 基本属性 ----------

        // 两条边向量（从顶点 a 出发）
        // 记住：几何计算永远用边向量，不是顶点坐标本身
        Vec3 edge_ab() const { return b - a; }
        Vec3 edge_ac() const { return c - a; }

        // 未归一化法向量（叉积，长度 = 2 * 面积）
        Vec3 rawNormal() const { return edge_ab().cross(edge_ac()); }

        // 面积 = 两边向量叉积长度的一半
        float area() const { return rawNormal().length() * 0.5f; }

        // 单位法向量
        Vec3 normal() const { return rawNormal().normalized(); }

        // 重心 = 三顶点坐标平均值
        Vec3 centroid() const { return (a + b + c) * (1.0f / 3.0f); }

        // ---------- 点与平面的关系 ----------

        // 点 P 到三角形所在平面的有符号距离
        // 正值 → P 在法向量正面一侧
        // 负值 → P 在法向量背面一侧
        float signedDistToPlane(const Vec3& P) const {
            Vec3 n = rawNormal();
            float len = n.length();
            if (len < EPS) return 0;
            return n.dot(P - a) / len;
        }

        // 判断点 P 是否在三角形所在平面内（给定误差 eps）
        bool pointOnPlane(const Vec3& P, float eps = 1e-5f) const {
            return std::abs(signedDistToPlane(P)) < eps;
        }

        // ---------- 点在三角形内的判断（3D）----------

        // 用重心坐标法判断点 P 是否在三角形内（含边界）
        // 前提：P 已经在三角形所在平面上（可先用 pointOnPlane 检查）
        // 返回：是否在三角形内，同时输出重心坐标 (u, v)，w = 1 - u - v
        bool inTriangle(const Vec3& P, float& u, float& v) const {
            Vec3 v0 = edge_ac();
            Vec3 v1 = edge_ab();
            Vec3 v2 = P - a;

            float dot00 = v0.dot(v0);
            float dot01 = v0.dot(v1);
            float dot02 = v0.dot(v2);
            float dot11 = v1.dot(v1);
            float dot12 = v1.dot(v2);

            float denom = dot00 * dot11 - dot01 * dot01;
            if (std::abs(denom) < EPS) return false;   // 退化三角形

            float inv = 1.0f / denom;
            u = (dot11 * dot02 - dot01 * dot12) * inv;
            v = (dot00 * dot12 - dot01 * dot02) * inv;

            return (u >= -EPS) && (v >= -EPS) && (u + v <= 1.0f + EPS);
        }

        // 不需要重心坐标时的简化版
        bool inTriangle(const Vec3& P) const {
            float u, v;
            return inTriangle(P, u, v);
        }

        // ---------- 光线与三角形求交（Möller–Trumbore 算法）----------
        // 用途：射线检测、碰撞检测，游戏/渲染引擎常用
        // 输入：光线起点 orig，方向 dir（不需要归一化）
        // 输出：t（交点在光线上的参数），u、v（重心坐标）
        // 返回：是否相交（t > 0 表示正方向相交）
        bool rayIntersect(const Vec3& orig, const Vec3& dir,
            float& t, float& u, float& v) const {
            Vec3 ab = edge_ab();
            Vec3 ac = edge_ac();
            Vec3 h = dir.cross(ac);
            float det = ab.dot(h);

            if (std::abs(det) < EPS) return false;   // 光线平行于三角形

            float inv_det = 1.0f / det;
            Vec3 s = orig - a;
            u = s.dot(h) * inv_det;
            if (u < 0 || u > 1) return false;

            Vec3 q = s.cross(ab);
            v = dir.dot(q) * inv_det;
            if (v < 0 || u + v > 1) return false;

            t = ac.dot(q) * inv_det;
            return t > EPS;   // t > 0 才是正方向上的交点
        }
    };


    // =============================================================================
    //  Part 5 : HalfEdge — 半边数据结构
    //
    //  半边结构是网格算法的核心，用于高效查询网格拓扑关系
    //
    //  核心概念：
    //    每条边拆成两条方向相反的半边（half-edge）
    //    三角形的三条半边按逆时针顺序排列
    //
    //  半边的三个关键指针：
    //    next  → 同一三角形内的下一条半边（逆时针）
    //    twin  → 对面的半边（相邻三角形的对应半边）
    //    vert  → 该半边的起点顶点索引
    //
    //  示意：
    //    三角形 [v0, v1, v2]
    //    半边 he0: v0→v1,  next=he1, twin=对面三角形的某条半边
    //    半边 he1: v1→v2,  next=he2
    //    半边 he2: v2→v0,  next=he0
    // =============================================================================
    struct Vertex
    {
        Vec3 pos;          // 顶点坐标
        int  outHE = -1;   // 任意一条从该顶点出发的半边索引

        Vertex() = default;
        explicit Vertex(const Vec3& p) : pos(p) {}
    };

    struct Face
    {
        int he = -1;    // 该面的任意一条半边索引

        Face() = default;
        explicit Face(int heIdx) : he(heIdx) {}

    };

    struct HalfEdge {
        int next = -1;   // 同面的下一条半边索引
        int twin = -1;   // 对面的半边索引（-1 表示边界边，没有相邻三角形）
        int vert = -1;   // 起点顶点索引
        int face = -1;   // 所属三角形索引（-1 表示边界虚半边）
    };

    // 半边网格（存储顶点 + 半边 + 面）
    struct HalfEdgeMesh {
        std::vector<Vertex>     verts;      // 顶点数组
        std::vector<Face>       faces;      // 面数组
        std::vector<HalfEdge>   halfedges;  // 半边数组

        // ---------- 构建 ----------

        // 从三角形列表构建半边结构
        // tris[i] = {v0, v1, v2}，顶点逆时针排列
        // Vec2 重载 — DelaunayMesh(2D网格)调用时使用，z 坐标补0
        void build(const std::vector<Vec2>& verts2d,
            const std::vector<std::array<int, 3>>& tris) {
            std::vector<Vec3> v3;
            v3.reserve(verts2d.size());
            for (const auto& v : verts2d)
                v3.push_back({ v.x, v.y, 0.0f });
            build(v3, tris);
        }

        void build(const std::vector<Vec3>& positions,
            const std::vector<std::array<int, 3>>& tris) {
            int nv = (int)positions.size();
            int nf = (int)tris.size();

            // 初始化顶点
            verts.clear();
            verts.reserve(nv);
            for (const auto& p : positions) {
                verts.emplace_back(p);
            }

            // 初始化面和半边
            faces.clear();
            faces.resize(nf);
            halfedges.clear();
            halfedges.resize(nf * 3);

            // 每个三角形生成 3 条半边
            for (int f = 0; f < nf; ++f) {
                faces[f].he = f * 3;
                for (int i = 0; i < 3; ++i) {
                    int he = f * 3 + i;
                    halfedges[he].vert = tris[f][i];
                    halfedges[he].next = f * 3 + (i + 1) % 3;
                    halfedges[he].face = f;
                    halfedges[he].twin = -1;
                }
            }

            // 建立顶点出边
            for (int he = 0; he < (int)halfedges.size(); ++he) {
                verts[halfedges[he].vert].outHE = he;
            }

            // 建立 twin 关系：边 (a→b) 的 twin 是 (b→a)
            std::map<std::pair<int, int>, int> edgeMap;
            for (int he = 0; he < (int)halfedges.size(); ++he) {
                int v0 = halfedges[he].vert;
                int v1 = halfedges[halfedges[he].next].vert;
                edgeMap[{v0, v1}] = he;
            }
            for (int he = 0; he < (int)halfedges.size(); ++he) {
                int v0 = halfedges[he].vert;
                int v1 = halfedges[halfedges[he].next].vert;
                auto it = edgeMap.find({ v1, v0 });
                if (it != edgeMap.end()) {
                    halfedges[he].twin = it->second;
                }
            }

            
        }

        // ---------- 拓扑查询 ----------

        // 顶点 v 的所有相邻顶点（一环邻域，逆时针顺序）
        std::vector<int> neighborVerts(int v) const {
            std::vector<int> result;
            int start = verts[v].outHE;
            if (start == -1) return result;

            int cur = start;
            bool isBoundary = false;

            // 正向遍历（沿 twin→next）
            do {
                result.push_back(halfedges[halfedges[cur].next].vert);
                int tw = halfedges[cur].twin;
                if (tw == -1) { isBoundary = true; break; }
                cur = halfedges[tw].next;
            } while (cur != start);

            // 边界顶点：反向遍历另一侧
            if (isBoundary) {
                cur = start;
                while (true) {
                    // prev = cur.next.next（三角形中 cur 的上一条半边）
                    int prev = halfedges[halfedges[cur].next].next;
                    int tw = halfedges[prev].twin;
                    if (tw == -1) {
                        // prev 是边界边，prev.vert 是最后一个相邻顶点
                        result.push_back(halfedges[prev].vert);
                        break;
                    }
                    // prev.vert 是这一侧的相邻顶点
                    result.push_back(halfedges[prev].vert);
                    cur = halfedges[tw].next;  // 反向前进
                }
            }

            return result;
        }

        // 顶点 v 的所有相邻面（一环面邻域，逆时针顺序）
        std::vector<int> neighborFaces(int v) const {
            std::vector<int> result;
            int start = verts[v].outHE;
            if (start == -1) return result;

            int cur = start;
            bool isBoundary = false;

            // 正向遍历
            do
            {
                if (halfedges[cur].face != -1) result.push_back(halfedges[cur].face);

                int tw = halfedges[cur].twin;
                if (tw == -1) {
                    isBoundary = true;
                    break;
                }
                cur = halfedges[tw].next;
            } while (cur != start);

            // 边界顶点：反向遍历另一侧
            if (isBoundary) {
                cur = start;
                while (true) {
                    int prev = halfedges[halfedges[cur].next].next;
                    int tw = halfedges[prev].twin;
                    if (tw == -1) break;  // 边界边，对应的面正向已收集，停止
                    cur = halfedges[tw].next;
                    if (halfedges[cur].face != -1)
                        result.push_back(halfedges[cur].face);
                }
            }

            return result;
        }

        // 判断顶点 v 是否在边界上
        bool isBoundaryVert(int v) const {
            int start = verts[v].outHE;
            if (start == -1) return true;
            int he = start;
            do
            {
                if (halfedges[he].twin == -1) return true;
                he = halfedges[halfedges[he].twin].next;
            } while (he != start);

            return false;
        }

        // 面 f 的三个顶点索引
        std::array<int, 3> faceVerts(int f) const {
            int he = faces[f].he;
            return {
                halfedges[he].vert,
                halfedges[halfedges[he].next].vert,
                halfedges[halfedges[halfedges[he].next].next].vert
            };
        }


    };




    // =============================================================================
    //  Part 6 : 点定位（Point Location）
    //
    //  给定一个三角网格和一个查询点 P，找出 P 在哪个三角形内
    //  这是 Delaunay 增量插入算法的核心操作
    // =============================================================================

    // 线性扫描点定位（简单版，O(n)）
    // 适合网格规模较小时使用
    // 返回三角形索引，-1 表示点在网格外
    inline int locatePoint_linear(
        const std::vector<Vec2>& verts,
        const std::vector<std::array<int, 3>>& tris,
        const Vec2& P)
    {
        for (int i = 0; i < (int)tris.size(); i++) {
            const Vec2& A = verts[tris[i][0]];
            const Vec2& B = verts[tris[i][1]];
            const Vec2& C = verts[tris[i][2]];
            if (inTriangle(A, B, C, P)) return i;
        }
        return -1;
    }

    // 跳步点定位（Walk，O(sqrt(n)) 平均）
    // 从一个起始三角形出发，每步朝 P 所在方向跳到相邻三角形
    // 适合 Delaunay 插入时使用（通常从上次插入的三角形出发）
    // 需要 HalfEdgeMesh 提供拓扑信息
    inline int locatePoint_walk(
        const std::vector<Vec2>& verts,
        const HalfEdgeMesh& mesh,
        const std::vector<std::array<int, 3>>& tris,
        const Vec2& P,
        int startFace = 0)
    {
        int f = startFace;
        int maxIter = (int)tris.size() + 10;   // 防止死循环

        for (int iter = 0; iter < maxIter; iter++) {
            auto [v0, v1, v2] = mesh.faceVerts(f);
            const Vec2& A = verts[v0];
            const Vec2& B = verts[v1];
            const Vec2& C = verts[v2];

            // 检查 P 是否在当前三角形内
            if (inTriangle(A, B, C, P)) return f;

            // 找到 P 在哪条边的外侧，跳到对面三角形
            int he = mesh.faces[f].he;
            bool jumped = false;
            for (int i = 0; i < 3; i++) {
                int va = mesh.halfedges[he].vert;
                int vb = mesh.halfedges[mesh.halfedges[he].next].vert;
                // P 在边 va→vb 的右侧，说明要往这个方向走
                if (cross2D(verts[va], verts[vb], P) < -EPS) {
                    int tw = mesh.halfedges[he].twin;
                    if (tw == -1) return -1;   // 边界，P 在网格外
                    f = mesh.halfedges[tw].face;
                    jumped = true;
                    break;
                }
                he = mesh.halfedges[he].next;
            }
            if (!jumped) return f;   // 没有需要跳的边，就在当前三角形内
        }
        return -1;   // 超出迭代次数
    }


    // =============================================================================
    //  Part 7 : Delaunay 三角剖分（框架）
    // =============================================================================

    // 占位：Delaunay 网格结构
    struct DelaunayMesh {
        std::vector<Vec2>            verts;    // 顶点
        std::vector<std::array<int, 3>> tris;  // 三角形（顶点索引，逆时针）
        HalfEdgeMesh                 topo;    // 拓扑结构

        // -----------------------------------------------------------------------
        // 构造超级三角形，包含所有输入点
        // 返回超级三角形的三个顶点（不加入 verts，由 build 统一管理）
        // -----------------------------------------------------------------------
        Triangle2D createSuperTriangle(const std::vector<Vec2>& points) {

            float x_min = points[0].x; float x_max = points[0].x; float y_min = points[0].y; float y_max = points[0].y;

            for (const auto& p : points) {
                if (p.x < x_min) x_min = p.x;
                if (p.x > x_max) x_max = p.x;
                if (p.y < y_min) y_min = p.y;
                if (p.y > y_max) y_max = p.y;
            }

            float x_mid = (x_min + x_max) * 0.5f;
            float y_mid = (y_min + y_max) * 0.5f;

            float dx = x_max - x_min;
            float dy = y_max - y_min;
            float delta = (dx + dy + 1.0f) * 10.0f;

            Vec2 v0 = { x_mid - 2.0f * delta, y_mid - delta };
            Vec2 v1 = { x_mid + 2.0f * delta, y_mid - delta };
            Vec2 v2 = { x_mid,                y_mid + 2.0f * delta };

            return Triangle2D(v0, v1, v2);
        }

        // -----------------------------------------------------------------------
        // 线性扫描点定位：返回包含点 p 的三角形索引，-1 表示未找到
        // -----------------------------------------------------------------------
        int locatePoint(const Vec2& p) {
            for (int i = 0; i < (int)tris.size(); ++i) {
                const Vec2& A = verts[tris[i][0]];
                const Vec2& B = verts[tris[i][1]];
                const Vec2& C = verts[tris[i][2]];
                if (Geo::inTriangle(A, B, C, p)) return i;
            }

            return -1;
        }

        // -----------------------------------------------------------------------
        // 找到外接圆包含点 p 的所有三角形索引（空腔 cavity）
        // 原理：Bowyer-Watson 算法的核心步骤
        // -----------------------------------------------------------------------
        std::vector<int> findCavity(const Vec2& p) {
            std::vector<int> cavity;
            for (int i = 0; i < (int)tris.size(); ++i) {
                const Vec2& A = verts[tris[i][0]];
                const Vec2& B = verts[tris[i][1]];
                const Vec2& C = verts[tris[i][2]];
                if (Geo::inCircumcircle(A, B, C, p))
                    cavity.push_back(i);
            }

            return cavity;
        }

        // -----------------------------------------------------------------------
        // 找到 cavity 的边界边
        // 边界边定义：只属于一个 cavity 三角形的边（另一侧不在 cavity 内）
        // 返回：每条边界边用 {顶点索引A, 顶点索引B} 表示
        // -----------------------------------------------------------------------
        std::vector<std::array<int, 2>> findBoundaryEdges(const std::vector<int>& cavity) {

            // 统计 cavity 内每条边出现的次数
            // 出现1次 = 边界边；出现2次 = 内部边（两侧都在 cavity 内，需删除）
            std::map<std::array<int, 2>, int> edgeCount;

            for (int fi : cavity) {
                for (int i = 0; i < 3; ++i) {
                    int va = tris[fi][i];
                    int vb = tris[fi][(i + 1) % 3];
                    // 统一边的方向（小索引在前），避免 {a,b} 和 {b,a} 被当作不同边
                    std::array<int, 2> edge = { std::min(va, vb), std::max(va, vb) };
                    edgeCount[edge]++;
                }
            }

            std::vector<std::array<int, 2>> boundary;
            for (auto& [edge, cnt] : edgeCount) {
                if (cnt == 1) boundary.push_back(edge);
            }

            return boundary;
        }

        // -----------------------------------------------------------------------
        // 插入新点 P：Bowyer-Watson 算法
        // 步骤：① 找空腔  ② 删除空腔三角形  ③ 用边界边和新点构成新三角形
        // -----------------------------------------------------------------------
        void insertPoint(const Vec2& P) {
            // Step 1: 找到所有外接圆包含 P 的三角形（空腔）
            std::vector<int> cavity = findCavity(P);

            if (cavity.empty()) return;

            // Step 2：找到空腔的边界边
            std::vector<std::array<int, 2>> boundary = findBoundaryEdges(cavity);

            // Step 3：删除空腔中的三角形（从后往前删，避免索引失效）
            std::sort(cavity.begin(), cavity.end(), std::greater<int>());
            for (int fi : cavity) {
                tris.erase(tris.begin() + fi);
            }

            // Step 4：加入新顶点
            int newV = (int)verts.size();
            verts.push_back(P);

            // Step 5：用每条边界边和新点构成新三角形（保证逆时针）
            for (auto& edge : boundary) {
                int va = edge[0], vb = edge[1];
                // 检查 va→vb→newV 是否逆时针，若不是则交换 va vb
                float cross = Geo::cross2D(verts[va], verts[vb], verts[newV]);
                if (cross < 0) std::swap(va, vb);
                tris.push_back({ va, vb, newV });
            }
        }

        // -----------------------------------------------------------------------
        // 从点集构建完整 Delaunay 三角剖分（Bowyer-Watson 算法）
        // -----------------------------------------------------------------------
        void build(const std::vector<Vec2>& points) {
            verts.clear();
            tris.clear();

            // Step 1：构造超级三角形，其顶点索引为 0 1 2
            Triangle2D st = createSuperTriangle(points);
            verts.push_back(st.a);
            verts.push_back(st.b);
            verts.push_back(st.c);
            tris.push_back({ 0, 1, 2 });

            // Step 2：逐点插入
            for (const auto& p : points) {
                insertPoint(p);
            }

            // Step 3：删除含超级三角形顶点（索引 0 1 2）的所有三角形
            tris.erase(
                std::remove_if(tris.begin(), tris.end(),
                    [](const std::array<int, 3>& t) {
                        return t[0] <= 2 || t[1] <= 2 || t[2] <= 2;
                    }),
                tris.end()
                        );

            // Step 4：删除超级三角形的顶点，修正后续顶点索引
            verts.erase(verts.begin(), verts.begin() + 3);
            for (auto& t : tris) {
                t[0] -= 3;
                t[1] -= 3;
                t[2] -= 3;
            }

            // Step 5：重建半边拓扑
            topo.build(verts, tris);
        }
    };

    // =============================================================================
    //  Part 8 : 可视化导出
    //
    //  两种格式：
    //  1. exportCSV  → 输出顶点和三角形到 .csv，配合 Python + matplotlib 绘制
    //  2. exportSVG  → 直接生成 SVG 文件，浏览器打开即可查看
    // =============================================================================

    // 导出为 CSV（配合 Python 绘制）
    // 生成两个文件：prefix_verts.csv 和 prefix_tris.csv
    inline void exportCSV(
        const std::vector<Vec2>& verts,
        const std::vector<std::array<int, 3>>& tris,
        const std::string& prefix
    ) {
        std::ofstream fv(prefix + "_verts.csv");
        fv << "x,y\n";
        for (const auto& v : verts) {
            fv << v.x << "," << v.y << "\n";
            fv.close();

            std::ofstream ft(prefix + "_tris.csv");
            ft << "i0, i1, i2\n";
            for (const auto& t : tris) {
                ft << t[0] << "," << t[1] << "," << t[2] << "\n";
            }
            ft.close();
        }
    }

    // 导出为 SVG（直接在浏览器打开即可查看）
    // showIndex=true 时显示顶点编号，svgSize 控制画布大小
    inline void exportSVG(
        const std::vector<Vec2>& verts,
        const std::vector<std::array<int, 3>>& tris,
        const std::string& filename,
        bool showIndex = false,
        int svgSize = 600
    ) {
        if (verts.empty()) return;

        // 计算包围盒
        float xmin = verts[0].x, xmax = verts[0].x;
        float ymin = verts[0].y, ymax = verts[0].y;
        for (const auto& v : verts) {
            if (v.x < xmin) xmin = v.x;
            if (v.x > xmax) xmax = v.x;
            if (v.y < ymin) ymin = v.y;
            if (v.y > ymax) ymax = v.y;
        }
        const float margin = svgSize * 0.08f;
        const float drawSize = svgSize - 2.0f * margin;
        float range = std::max(xmax - xmin, ymax - ymin);
        if (range < 1e-9f) range = 1.0f;
        // 世界坐标 → SVG 坐标（y 轴翻转）
        auto sx = [&](float x) { return margin + (x - xmin) / range * drawSize; };
        auto sy = [&](float y) { return margin + (ymax - y) / range * drawSize; };

        std::ofstream f(filename);
        f << std::fixed;
        f.precision(2);

        // 使用单独的语句输出，避免字符串内嵌引号问题
        const std::string Q = "\"";
        f << "<svg xmlns=" << Q << "http://www.w3.org/2000/svg" << Q
            << " width=" << Q << svgSize << Q
            << " height=" << Q << svgSize << Q << ">\n";
        f << "<rect width=" << Q << "100%" << Q
            << " height=" << Q << "100%" << Q
            << " fill=" << Q << "#1e272e" << Q << "/>\n";

        // 三角形边
        f << "<g stroke=" << Q << "#74b9ff" << Q
            << " stroke-width=" << Q << "0.8" << Q
            << " fill=" << Q << "none" << Q
            << " opacity=" << Q << "0.7" << Q << ">\n";
        for (const auto& t : tris) {
            float ax = sx(verts[t[0]].x), ay = sy(verts[t[0]].y);
            float bx = sx(verts[t[1]].x), by = sy(verts[t[1]].y);
            float cx = sx(verts[t[2]].x), cy = sy(verts[t[2]].y);
            f << "<polygon points=" << Q
                << ax << "," << ay << " "
                << bx << "," << by << " "
                << cx << "," << cy << Q << "/>\n";
        }
        f << "</g>\n";

        // 顶点（黄色圆点）
        f << "<g fill=" << Q << "#fdcb6e" << Q
            << " stroke=" << Q << "#e17055" << Q
            << " stroke-width=" << Q << "0.5" << Q << ">\n";
        for (const auto& v : verts) {
            f << "<circle cx=" << Q << sx(v.x) << Q
                << " cy=" << Q << sy(v.y) << Q
                << " r=" << Q << "3" << Q << "/>\n";
        }
        f << "</g>\n";

        // 顶点编号（可选）
        if (showIndex) {
            f << "<g font-size=" << Q << "9" << Q
                << " fill=" << Q << "#dfe6e9" << Q
                << " font-family=" << Q << "monospace" << Q << ">\n";
            for (int i = 0; i < (int)verts.size(); ++i) {
                f << "<text x=" << Q << (sx(verts[i].x) + 4) << Q
                    << " y=" << Q << (sy(verts[i].y) - 4) << Q
                    << ">" << i << "</text>\n";
            }
            f << "</g>\n";
        }

        // 统计信息
        f << "<text x=" << Q << "10" << Q
            << " y=" << Q << "20" << Q
            << " font-size=" << Q << "12" << Q
            << " fill=" << Q << "#95a5a6" << Q
            << " font-family=" << Q << "Arial" << Q << ">"
            << "verts=" << verts.size() << "  tris=" << tris.size()
            << "</text>\n";

        f << "</svg>\n";
        f.close();
    }

    // DelaunayMesh 的便捷导出方法
    inline void exportMeshCSV(const DelaunayMesh& mesh, const std::string& prefix) {
        exportCSV(mesh.verts, mesh.tris, prefix);
    }

    inline void exportMeshSVG(const DelaunayMesh& mesh, const std::string& filename,
        bool showIndex = false, int svgSize = 600) {
        exportSVG(mesh.verts, mesh.tris, filename, showIndex, svgSize);
    }

}   // namespace Geo