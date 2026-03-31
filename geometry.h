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
//  2026-03  Part 4.5: 网格拓扑检验 (isManifold)
//  2026-03  Part 4.6: 几何退化检测 (Degeneracy Check)
//  2026-03  Part 4.7: 自相交检测 (Self-Intersection)
// =============================================================================


// =============================================================================
//  Part 0 : 常量与工具函数
// =============================================================================

namespace Geo {

    // 浮点类型别名——方便切换 float/double
    // Delaunay 等算法对数值精度敏感，推荐使用 double
    using real = double;

    constexpr real EPS = 1e-9;	 // 浮点判零阈值，几何判断用
    constexpr real PI = 3.14159265358979323846;

    // 浮点判零
    inline bool isZero(real x) { return std::abs(x) < EPS; }

    // 浮点判等
    inline bool equal(real a, real b) { return std::abs(a - b) < EPS; }

    // 符号函数：返回 -1 / 0 / 1
    inline int sign(real x) {
        if (x > EPS) return  1;
        if (x < -EPS) return -1;
        return 0;
    }


    // =============================================================================
    //  Part 1 : Vec2 — 二维向量
    // =============================================================================

    struct Vec2 {
        real x, y;

        // ---------- 构造 ----------
        Vec2() : x(0), y(0) {}
        Vec2(real x, real y) : x(x), y(y) {}

        // ---------- 基本运算 ----------
        Vec2 operator+(const Vec2& v) const { return { x + v.x, y + v.y }; }
        Vec2 operator-(const Vec2& v) const { return { x - v.x, y - v.y }; }
        Vec2 operator*(real f)       const { return { x * f,   y * f }; }
        Vec2 operator/(real f)       const { assert(std::abs(f) > EPS && "Vec2: division by zero"); return { x / f,   y / f }; }

        Vec2& operator+=(const Vec2& v) { x += v.x; y += v.y; return *this; }
        Vec2& operator-=(const Vec2& v) { x -= v.x; y -= v.y; return *this; }
        Vec2& operator*=(real f) { x *= f;   y *= f;   return *this; }

        bool operator==(const Vec2& v) const { return equal(x, v.x) && equal(y, v.y); }
        bool operator!=(const Vec2& v) const { return !(*this == v); }

        // ---------- 向量运算 ----------

        // 点积
        real dot(const Vec2& v) const { return x * v.x + y * v.y; }

        // 2D叉积（标量）
        // 正值 → v 在 *this 的逆时针方向
        // 负值 → v 在 *this 的顺时针方向
        // 零   → 平行或共线
        real cross(const Vec2& v) const { return x * v.y - y * v.x; }

        // 向量长度
        real length()  const { return std::sqrt(x * x + y * y); }
        real length2() const { return x * x + y * y; }   // 长度的平方，避免开方更快

        // 单位向量（返回新向量，不修改自身）
        // 零向量时返回自身，不崩溃
        Vec2 normalized() const {
            real len = length();
            if (len < EPS) return *this;
            return *this / len;
        }

        // 逆时针旋转 90 度的垂直向量
        Vec2 perp() const { return { -y, x }; }

        // ---------- 工具 ----------
        bool isZero() const { return length2() < EPS; }
    };

    // 支持 scalar * Vec2 写法
    inline Vec2 operator*(real f, const Vec2& v) { return v * f; }

    // 两点距离
    inline real dist(const Vec2& a, const Vec2& b) { return (a - b).length(); }
    inline real dist2(const Vec2& a, const Vec2& b) { return (a - b).length2(); }


    // =============================================================================
    //  Part 2 : Vec3 — 三维向量
    // =============================================================================

    struct Vec3 {
        real x, y, z;

        // ---------- 构造 ----------
        Vec3() : x(0), y(0), z(0) {}
        Vec3(real x, real y, real z) : x(x), y(y), z(z) {}

        // ---------- 基本运算 ----------
        Vec3 operator+(const Vec3& v) const { return { x + v.x, y + v.y, z + v.z }; }
        Vec3 operator-(const Vec3& v) const { return { x - v.x, y - v.y, z - v.z }; }
        Vec3 operator*(real f)       const { return { x * f,   y * f,   z * f }; }
        Vec3 operator/(real f)       const { assert(std::abs(f) > EPS && "Vec3: division by zero"); return { x / f,   y / f,   z / f }; }

        Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
        Vec3& operator-=(const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
        Vec3& operator*=(real f) { x *= f;   y *= f;   z *= f;   return *this; }

        bool operator==(const Vec3& v) const { return equal(x, v.x) && equal(y, v.y) && equal(z, v.z); }
        bool operator!=(const Vec3& v) const { return !(*this == v); }

        // ---------- 向量运算 ----------

        // 点积
        real dot(const Vec3& v) const {
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
        real length()  const { return std::sqrt(x * x + y * y + z * z); }
        real length2() const { return x * x + y * y + z * z; }

        // 单位向量（零向量时返回自身）
        Vec3 normalized() const {
            real len = length();
            if (len < EPS) return *this;
            return *this / len;
        }

        // ---------- 工具 ----------
        bool isZero() const { return length2() < EPS; }

        // 投影到另一个向量上的标量分量
        real projectOnto(const Vec3& axis) const {
            return this->dot(axis.normalized());
        }
    };

    inline Vec3 operator*(real f, const Vec3& v) { return v * f; }

    inline real dist(const Vec3& a, const Vec3& b) { return (a - b).length(); }
    inline real dist2(const Vec3& a, const Vec3& b) { return (a - b).length2(); }

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
    inline real cross2D(const Vec2& O, const Vec2& A, const Vec2& B) {
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
        real d1 = cross2D(C, D, A);
        real d2 = cross2D(C, D, B);
        real d3 = cross2D(A, B, C);
        real d4 = cross2D(A, B, D);

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
    inline real triangleSignedArea(const Vec2& A, const Vec2& B, const Vec2& C) {
        return cross2D(A, B, C) * 0.5;
    }

    // 三角形面积（始终为正）
    inline real triangleArea(const Vec2& A, const Vec2& B, const Vec2& C) {
        return std::abs(triangleSignedArea(A, B, C));
    }

    // 判断点 P 是否在三角形 ABC 内（包括边上）
    // 原理：P 在 ABC 内 ⟺ P 在 AB、BC、CA 三条有向边的同侧
    // 注意：无论 ABC 顶点是顺时针还是逆时针，has_neg && has_pos 的判断都正确
    inline bool inTriangle(const Vec2& A, const Vec2& B, const Vec2& C, const Vec2& P) {
        real d1 = cross2D(A, B, P);
        real d2 = cross2D(B, C, P);
        real d3 = cross2D(C, A, P);
        bool has_neg = (d1 < -EPS) || (d2 < -EPS) || (d3 < -EPS);
        bool has_pos = (d1 > EPS) || (d2 > EPS) || (d3 > EPS);
        return !(has_neg && has_pos);
    }

    // 三角形外接圆圆心
    // 原理：外接圆圆心是三边垂直平分线的交点
    // 用途：Delaunay 三角剖分的核心计算
    inline Vec2 circumcenter(const Vec2& A, const Vec2& B, const Vec2& C) {
        real ax = B.x - A.x, ay = B.y - A.y;
        real bx = C.x - A.x, by = C.y - A.y;
        real D = 2.0 * (ax * by - ay * bx);
        // D 接近零说明三点共线，外接圆不存在
        assert(!isZero(D) && "circumcenter: three points are collinear");
        real ux = (by * (ax * ax + ay * ay) - ay * (bx * bx + by * by)) / D;
        real uy = (ax * (bx * bx + by * by) - bx * (ax * ax + ay * ay)) / D;
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
        real area2 = cross2D(A, B, C);
        if (isZero(area2)) return false;

        // 如果顶点是顺时针排列，交换 B C 使其变为逆时针
        // 保证 det 的符号含义一致
        const Vec2* a = &A;
        const Vec2* b = &B;
        const Vec2* c = &C;
        if (area2 < 0) std::swap(b, c);

        // 平移：以 P 为原点，A B C 坐标相应减去 P
        real ax = a->x - P.x, ay = a->y - P.y;
        real bx = b->x - P.x, by = b->y - P.y;
        real cx = c->x - P.x, cy = c->y - P.y;

        // 各点到 P 的距离平方（提升到抛物面的 z 坐标）
        real ar2 = ax * ax + ay * ay;
        real br2 = bx * bx + by * by;
        real cr2 = cx * cx + cy * cy;

        // 3×3 行列式展开（按第三列余子式）
        real det = ax * (by * cr2 - cy * br2)
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
        Triangle2D() : a(), b(), c() {}
        Triangle2D(const Vec2& a, const Vec2& b, const Vec2& c) : a(a), b(b), c(c) {}

        // ---------- 基本属性 ----------

        // 边向量（从顶点 a 出发）
        Vec2 edge_ab() const { return b - a; }
        Vec2 edge_ac() const { return c - a; }

        // 有符号面积（逆时针为正）
        real signedArea() const { return triangleSignedArea(a, b, c); }

        // 面积（始终为正）
        real area() const { return triangleArea(a, b, c); }

        // 重心
        Vec2 centroid() const { return (a + b + c) * (1.0 / 3.0); }

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
    inline real polygonSignedArea(const std::vector<Vec2>& poly) {
        real area = 0;
        int n = (int)poly.size();
        for (int i = 0; i < n; i++) {
            const Vec2& cur = poly[i];
            const Vec2& next = poly[(i + 1) % n];
            area += cur.cross(next);
        }
        return area * 0.5;
    }

    inline real polygonArea(const std::vector<Vec2>& poly) {
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
            real yi = poly[i].y, yj = poly[j].y;
            real xi = poly[i].x, xj = poly[j].x;
            if ((yi > P.y) != (yj > P.y)) {
                real xIntersect = (xj - xi) * (P.y - yi) / (yj - yi) + xi;
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
            real c = cross2D(p0, a, b);
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

    // 辅助结构：三维轴向包围盒 (Axis-Aligned Bounding Box)
    struct AABB3D
    {
        Vec3 min_pt, max_pt;

        AABB3D() = default;
        AABB3D(const Vec3& a, const Vec3& b, const Vec3& c) {
            min_pt.x = std::min({ a.x, b.x, c.x });
            min_pt.y = std::min({ a.y, b.y, c.y });
            min_pt.z = std::min({ a.z, b.z, c.z });

            max_pt.x = std::max({ a.x, b.x, c.x });
            max_pt.y = std::max({ a.y, b.y, c.y });
            max_pt.z = std::max({ a.z, b.z, c.z });
        }

        // 检查两个包围盒是否发生空间重叠
        bool intersect(const AABB3D& other) const {
            if (max_pt.x < other.min_pt.x || min_pt.x > other.max_pt.x) return false;
            if (max_pt.y < other.min_pt.y || min_pt.y > other.max_pt.y) return false;
            if (max_pt.z < other.min_pt.z || min_pt.z > other.max_pt.z) return false;
            return true;
        }
    };

    struct Triangle3D {
        Vec3 a, b, c;   // 三个顶点，建议逆时针朝向正面

        // ---------- 构造 ----------
        Triangle3D() : a(), b(), c() {}
        Triangle3D(const Vec3& a, const Vec3& b, const Vec3& c) : a(a), b(b), c(c) {}

        // ---------- 基本属性 ----------

        // 两条边向量（从顶点 a 出发）
        // 记住：几何计算永远用边向量，不是顶点坐标本身
        Vec3 edge_ab() const { return b - a; }
        Vec3 edge_ac() const { return c - a; }

        // 未归一化法向量（叉积，长度 = 2 * 面积）
        Vec3 rawNormal() const { return edge_ab().cross(edge_ac()); }

        // 面积 = 两边向量叉积长度的一半
        real area() const { return rawNormal().length() * 0.5; }

        // 单位法向量
        Vec3 normal() const { return rawNormal().normalized(); }

        // 重心 = 三顶点坐标平均值
        Vec3 centroid() const { return (a + b + c) * (1.0 / 3.0); }

        // ---------- 点与平面的关系 ----------

        // 点 P 到三角形所在平面的有符号距离
        // 正值 → P 在法向量正面一侧
        // 负值 → P 在法向量背面一侧
        real signedDistToPlane(const Vec3& P) const {
            Vec3 n = rawNormal();
            real len = n.length();
            if (len < EPS) return 0;
            return n.dot(P - a) / len;
        }

        // 判断点 P 是否在三角形所在平面内（给定误差 eps）
        bool pointOnPlane(const Vec3& P, real eps = 1e-5) const {
            return std::abs(signedDistToPlane(P)) < eps;
        }

        // ---------- 点在三角形内的判断（3D）----------

        // 用重心坐标法判断点 P 是否在三角形内（含边界）
        // 前提：P 已经在三角形所在平面上（可先用 pointOnPlane 检查）
        // 返回：是否在三角形内，同时输出重心坐标 (u, v)，w = 1 - u - v
        bool inTriangle(const Vec3& P, real& u, real& v) const {
            Vec3 v0 = edge_ac();
            Vec3 v1 = edge_ab();
            Vec3 v2 = P - a;

            real dot00 = v0.dot(v0);
            real dot01 = v0.dot(v1);
            real dot02 = v0.dot(v2);
            real dot11 = v1.dot(v1);
            real dot12 = v1.dot(v2);

            real denom = dot00 * dot11 - dot01 * dot01;
            if (std::abs(denom) < EPS) return false;   // 退化三角形

            real inv = 1.0 / denom;
            u = (dot11 * dot02 - dot01 * dot12) * inv;
            v = (dot00 * dot12 - dot01 * dot02) * inv;

            return (u >= -EPS) && (v >= -EPS) && (u + v <= 1.0 + EPS);
        }

        // 不需要重心坐标时的简化版
        bool inTriangle(const Vec3& P) const {
            real u, v;
            return inTriangle(P, u, v);
        }

        // ---------- 光线与三角形求交（Möller–Trumbore 算法）----------
        // 用途：射线检测、碰撞检测，游戏/渲染引擎常用
        // 输入：光线起点 orig，方向 dir（不需要归一化）
        // 输出：t（交点在光线上的参数），u、v（重心坐标）
        // 返回：是否相交（t > 0 表示正方向相交）
        bool rayIntersect(const Vec3& orig, const Vec3& dir,
            real& t, real& u, real& v) const {
            Vec3 ab = edge_ab();
            Vec3 ac = edge_ac();
            Vec3 h = dir.cross(ac);
            real det = ab.dot(h);

            if (std::abs(det) < EPS) return false;   // 光线平行于三角形

            real inv_det = 1.0 / det;
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
    //  Part 4.5 : 网格拓扑检验 (Mesh Validation)
    //
    //  在构建 HalfEdgeMesh 等高级拓扑结构前，必须确保网格是流形 (Manifold)。
    //  包含：检测非流形边、检测非流形顶点。
    // =============================================================================

    // 检查Polygon Soup 是否是合法的流形网格
    inline bool isManifold(const std::vector<std::array<int, 3>>& tris) {
        if (tris.empty()) return true;

        // Step 0: 获取最大顶点索引，用于动态分配邻接表内存
        int max_v = -1;
        for (const auto& t : tris) {
            max_v = std::max({ max_v, t[0], t[1], t[2] });
        }
        if (max_v == -1) return true;

        // Step 1: 检测非流形边 (Non-manifold Edges)
        // 规则：任何一条无向边，最多只能被 2 个面共享
        std::vector<std::array<int, 2>> edges;
        edges.reserve(tris.size() * 3);

        for (const auto& t : tris) {
            for (int i = 0; i < 3; ++i) {
                int u = t[i];
                int v = t[(i + 1) % 3];
                edges.push_back({ std::min(u, v), std::max(u, v) });
            }
        }

        // 排序，将相同的边挨在一起
        std::sort(edges.begin(), edges.end());

        int n_edges = (int)edges.size();
        for (int i = 0; i < n_edges;) {
            int count = 1;
            while (i + count < n_edges && edges[i][0] == edges[i + count][0]
                                       && edges[i][1] == edges[i + count][1])
            {
                count++;
            }

            if (count > 2) {
                return false;
            }
            i += count;
        }

        // Step 2: 检测非流形顶点 (Non-manifold Vertices)
        // 规则：围绕任何一个顶点的所有面，必须能通过共享边连成一个单一的连通块
        std::vector<std::vector<int>> vert_to_faces(max_v + 1);
        for (int i = 0; i < (int)tris.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                vert_to_faces[tris[i][j]].push_back(i);
            }
        }

        for (int v = 0; v <= max_v; ++v) {
            const auto& faces = vert_to_faces[v];
            int n_faces = (int)faces.size();

            if (n_faces <= 1) continue;

            std::vector<bool> visited(n_faces, false);
            std::vector<int> queue;
            queue.reserve(n_faces);

            visited[0] = true;
            queue.push_back(0);
            int head = 0;
            int visited_count = 1;

            while (head < (int)queue.size()) {
                int curr_face_idx = queue[head++];
                int curr_face = faces[curr_face_idx];

                // 扫描周围其他未访问的面，看是否与 curr_face 相邻
                for (int i = 0; i < n_faces; ++i) {
                    if (visited[i]) continue;
                    int neighbor_face = faces[i];

                    // 检查判定：curr_face 和 neighbor_face 是否共享了除 v 以外的另一个顶点？
                    // 如果共享了，说明它们在 v 旁边紧紧挨着（共享了一条包含 v 的边）
                    int shared_other_v = 0;
                    for (int j = 0; j < 3; ++j) {
                        int v1 = tris[curr_face][j];
                        if (v1 == v) continue;

                        for (int k = 0; k < 3; ++k) {
                            int v2 = tris[neighbor_face][k];
                            if (v2 == v) continue;

                            if (v1 == v2) {
                                shared_other_v++;
                            }
                        }
                    }

                    // 找到了相连的面，标记为访问并推入队列
                    if (shared_other_v >= 1) {
                        visited[i] = true;
                        queue.push_back(i);
                        visited_count++;
                    }
                }
            }
            // BFS 扩散结束。如果最终访问的面数量小于这个顶点连接的总面数，
            // 说明在顶点 v 处存在互不相连的"沙漏"结构面簇。
            if (visited_count < n_faces) {
                return false; // 触发非流形顶点
            }
        }

        return true;
    }

    // =============================================================================
    //  Part 4.6 : 几何退化检测 (Geometric Degeneracy Validation)
    // =============================================================================

    // 找出所有退化的三角形（面积极小，或存在极大/极小内角）
    // 输入：3D 顶点数组 verts，三角形索引数组 tris
    // 返回：退化三角形的索引列表
    inline std::vector<int> findDegenerateTriangles(
        const std::vector<Vec3>& verts,
        const std::vector<std::array<int, 3>>& tris,
        real minArea = 1e-7,       // 最小允许面积
        real minAngleDeg = 1.0,    // 最小允许内角（度）—— 拦截“针状 (Needle)”
        real maxAngleDeg = 179.0   // 最大允许内角（度）—— 拦截“帽状 (Cap)”
    ) {
        std::vector<int> degenerate_indices;

        real minCos = std::cos(maxAngleDeg * PI / 180.0);
        real maxCos = std::cos(minAngleDeg * PI / 180.0);

        for (int i = 0; i < (int)tris.size(); ++i) {
            const Vec3& A = verts[tris[i][0]];
            const Vec3& B = verts[tris[i][1]];
            const Vec3& C = verts[tris[i][2]];

            Vec3 ab = B - A;
            Vec3 ac = C - A;
            Vec3 bc = C - B;

            // 拦截 1：面积检测 (Area Check)
            // 用叉积求面积：0.5 * ||AB x AC||
            real area = ab.cross(ac).length() * 0.5;
            if (area < minArea) {
                degenerate_indices.push_back(i);
                continue;   // 已经是退化面，直接跳过后续角度计算
            }

            // 拦截 2：角度检测 (Angle Check / Aspect Ratio Check)
            real len_ab = ab.length();
            real len_ac = ac.length();
            real len_bc = bc.length();

            // 防御性编程：理论上面积检测已经拦截了极短边，但防止浮点截断导致除以 0
            if (len_ab < EPS || len_ac < EPS || len_bc < EPS) {
                degenerate_indices.push_back(i);
                continue;
            }

            // 算出从每个顶点出发的单位方向向量
            Vec3 dir_ab = ab / len_ab;
            Vec3 dir_ac = ac / len_ac;
            Vec3 dir_ba = ab * -1.0 / len_ab;
            Vec3 dir_bc = bc / len_bc;
            Vec3 dir_ca = ac * -1.0 / len_ac;
            Vec3 dir_cb = bc * -1.0 / len_bc;

            // 用向量点积求三个内角的余弦值 (Dot Product)
            // cos(theta) = Dir1 · Dir2
            real cosA = dir_ab.dot(dir_ac);
            real cosB = dir_ba.dot(dir_bc);
            real cosC = dir_ca.dot(dir_cb);

            // 判决逻辑：
            // 如果 cos 值 > maxCos (接近 1)，说明夹角极其尖锐 (小于 minAngleDeg)
            // 如果 cos 值 < minCos (接近 -1)，说明夹角极其平缓 (大于 maxAngleDeg)
            if (cosA > maxCos || cosA < minCos ||
                cosB > maxCos || cosB < minCos ||
                cosC > maxCos || cosC < minCos) {
                degenerate_indices.push_back(i);
            }
        }

        return degenerate_indices;
    }

    // =============================================================================
    //  Part 4.7 : 自相交检测 (Self-Intersection Detection)
    // =============================================================================

    // 全局自相交检测函数
    // 输入：3D 顶点数组 verts，三角形索引数组 tris
    // 返回：所有发生互相穿透的三角形索引对 (Tri_A, Tri_B)
    inline std::vector<std::pair<int, int>> findSelfIntersections(
        const std::vector<Vec3>& verts,
        const std::vector<std::array<int, 3>>& tris
    ) {
        std::vector<std::pair<int, int>> intersecting_pairs;
        int n_tris = (int)tris.size();

        // 准备阶段：预先计算所有面片的 AABB 和 Triangle3D 对象 (空间换时间)
        // 避免在 O(N^2) 的双重循环中反复进行内存分配和重复计算
        std::vector<AABB3D> aabbs;
        std::vector<Triangle3D> tri_objs;
        aabbs.reserve(n_tris);
        tri_objs.reserve(n_tris);

        for (int i = 0; i < n_tris; ++i) {
            const Vec3& A = verts[tris[i][0]];
            const Vec3& B = verts[tris[i][1]];
            const Vec3& C = verts[tris[i][2]];
            aabbs.emplace_back(A, B, C);
            tri_objs.emplace_back(A, B, C);
        }

        // 核心检测阶段：双重循环遍历所有可能的三角形对
        for (int i = 0; i < n_tris; ++i) {
            for (int j = i + 1; j < n_tris; ++j) {
                // 拦截 1：拓扑邻接过滤 (Adjacency Filter)
                // 如果两个三角形共享了 1 个、2 个或 3 个顶点，说明它们在拓扑上连在一起。
                // 它们在边界处的物理接触是合法的，不属于“自相交穿透”。
                int shared_verts = 0;
                for (int vi = 0; vi < 3; ++vi) {
                    for (int vj = 0; vj < 3; ++vj) {
                        if (tris[i][vi] == tris[j][vj]) {
                            shared_verts++;
                        }
                    }
                }
                if (shared_verts > 0) continue;

                // 拦截 2：空间粗筛 (AABB Pre-rejection)
                if (!aabbs[i].intersect(aabbs[j])) continue;

                // 拦截 3：精确几何求交 (Edge-Piercing Test)
                const Triangle3D& T1 = tri_objs[i];
                const Triangle3D& T2 = tri_objs[j];
                bool is_intersect = false;

                // 辅助 Lambda：测试线段 p0->p1 是否刺穿了三角形 T
                auto segmentPiercesTri = [](const Vec3& p0, const Vec3& p1, const Triangle3D& T) {
                    real t, u, v;
                    Vec3 dir = p1 - p0;

                    // 调用在 Part 4 中已经实现的 Möller–Trumbore 射线求交算法
                    if (T.rayIntersect(p0, dir, t, u, v)) {
                        // 射线求交返回的是 t，线段的范围是 t 属于 (0, 1)
                        // 使用 EPS 缩紧边界，防止极端贴合导致的误判
                        if (t > EPS && t < 1.0 - EPS) {
                            return true;
                        }
                    }
                    return false;
                };

                // 测试 T1 的三条边是否刺穿 T2
                if (segmentPiercesTri(T1.a, T1.b, T2) ||
                    segmentPiercesTri(T1.b, T1.c, T2) ||
                    segmentPiercesTri(T1.c, T1.a, T2)) {
                    is_intersect = true;
                }
                // 测试 T2 的三条边是否刺穿 T1
                else if (segmentPiercesTri(T2.a, T2.b, T1) ||
                    segmentPiercesTri(T2.b, T2.c, T1) ||
                    segmentPiercesTri(T2.c, T2.a, T1)) {
                    is_intersect = true;
                }

                if (is_intersect) {
                    intersecting_pairs.push_back({ i, j });
                }
            }
        }

        return intersecting_pairs;
    }


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
                v3.push_back({ v.x, v.y, 0.0 });
            build(v3, tris);
        }

        void build(const std::vector<Vec3>& positions,
            const std::vector<std::array<int, 3>>& tris) {

            // 防御性拦截：确保输入网格拓扑合法，避免半边结构指针错乱崩溃
            assert(Geo::isManifold(tris) && "HalfEdgeMesh requires manifold topology!");

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
            // 优化：使用 std::vector 替代 std::map 以提升 CPU 缓存命中率和百万级面片加载速度
            struct DirectedEdge {
                int v0, v1;
                int he_idx;
                bool operator<(const DirectedEdge& other) const {
                    if (v0 != other.v0) return v0 < other.v0;
                    return v1 < other.v1;
                }
            };

            std::vector<DirectedEdge> edgeList;
            edgeList.reserve(halfedges.size());
            for (int he = 0; he < (int)halfedges.size(); ++he) {
                int v0 = halfedges[he].vert;
                int v1 = halfedges[halfedges[he].next].vert;
                edgeList.push_back({ v0, v1, he });
            }

            // 排序，使起点和终点相同的边在内存中连续排列
            std::sort(edgeList.begin(), edgeList.end());

            for (int he = 0; he < (int)halfedges.size(); ++he) {
                int v0 = halfedges[he].vert;
                int v1 = halfedges[halfedges[he].next].vert;

                // 使用二分查找寻找对向边 (v1 → v0)
                DirectedEdge target = { v1, v0, -1 };
                auto it = std::lower_bound(edgeList.begin(), edgeList.end(), target);

                if (it != edgeList.end() && it->v0 == v1 && it->v1 == v0) {
                    halfedges[he].twin = it->he_idx;
                }
            }
        }



        // ---------- 拓扑查询 ----------

        // 辅助：找到顶点 v 最右侧的边界出边
        // 如果 v 是内部顶点（没有边界边），返回 -1
        int findRightmostHE(int v) const {
            int start = verts[v].outHE;
            if (start == -1) return -1;
            int cur = start;
            do {
                int tw = halfedges[cur].twin;
                if (tw == -1) return cur;
                cur = halfedges[tw].next;
            } while (cur != start);
            return -1;  // 绕了一圈，内部顶点
        }

        // 顶点 v 的所有相邻顶点（一环邻域）
        // 统一保证返回的顶点顺序为严格的逆时针 (CCW)
        std::vector<int> neighborVerts(int v) const {
            std::vector<int> result;
            int start = verts[v].outHE;
            if (start == -1) return result;

            int rightmost = findRightmostHE(v);

            if (rightmost == -1) {
                // 内部顶点：统一改为逆时针 (CCW) 遍历
                int cur = start;
                do {
                    result.push_back(halfedges[halfedges[cur].next].vert);
                    int prev = halfedges[halfedges[cur].next].next;
                    cur = halfedges[prev].twin; // 逆时针步进
                } while (cur != start);
            }
            else {
                // 边界顶点：从最右侧边界出边逆时针 (CCW) 遍历
                int cur = rightmost;
                result.push_back(halfedges[halfedges[cur].next].vert);
                while (true) {
                    int prev = halfedges[halfedges[cur].next].next;
                    result.push_back(halfedges[prev].vert);

                    int tw = halfedges[prev].twin;
                    if (tw == -1) break; // 碰到最左侧边界
                    cur = tw;
                }
            }
            return result;
        }

        // 顶点 v 的所有相邻面（一环面邻域）
        // 统一保证返回的面顺序为严格的逆时针 (CCW)
        std::vector<int> neighborFaces(int v) const {
            std::vector<int> result;
            int start = verts[v].outHE;
            if (start == -1) return result;

            int rightmost = findRightmostHE(v);

            if (rightmost == -1) {
                // 内部顶点：统一改为逆时针 (CCW) 遍历
                int cur = start;
                do {
                    result.push_back(halfedges[cur].face);
                    int prev = halfedges[halfedges[cur].next].next;
                    cur = halfedges[prev].twin;
                } while (cur != start);
            }
            else {
                // 边界顶点：从最右侧边界出边逆时针 (CCW) 遍历
                int cur = rightmost;
                while (true) {
                    result.push_back(halfedges[cur].face);

                    int prev = halfedges[halfedges[cur].next].next;
                    int tw = halfedges[prev].twin;
                    if (tw == -1) break;
                    cur = tw;
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

        // 边翻转
        void flipEdge(int he) {
            int e = he;                                 // D→B
            int e_n = halfedges[e].next;                // B→A
            int e_p = halfedges[e_n].next;              // A→D

            int et = halfedges[e].twin;                 // B→D
            int et_n = halfedges[et].next;              // D→C
            int et_p = halfedges[et_n].next;            // C→B

            int f1 = halfedges[e].face;                 // ADB 改成 ACB
            int f2 = halfedges[et].face;                // BDC 改成 ADC

            int A = halfedges[e_p].vert;
            int C = halfedges[et_p].vert;

            // 修改翻转边起点
            halfedges[e].vert = A;
            halfedges[et].vert = C;

            // 重连 next
            halfedges[e].next = et_p;
            halfedges[et_p].next = e_n;
            halfedges[e_n].next = e;

            halfedges[et].next = e_p;
            halfedges[e_p].next = et_n;
            halfedges[et_n].next = et;

            // 修改 face 归属
            halfedges[e].face = f1;
            halfedges[et_p].face = f1;
            halfedges[e_n].face = f1;

            halfedges[e_p].face = f2;
            halfedges[et_n].face = f2;
            halfedges[et].face = f2;

            // 更新 Face 的代表半边
            faces[f1].he = e;
            faces[f2].he = e_p;
        }

        // 边分裂
        // 注意：只能分裂内部边（twin != -1），边界边不可分裂
        int splitEdge(int he) {
            int e = he;

            // 边界检查：边界边没有对面三角形，无法分裂
            if (halfedges[e].twin == -1) {
                assert(false && "splitEdge: cannot split a boundary edge");
                return -1;
            }

            int e_n = halfedges[e].next;
            int e_p = halfedges[e_n].next;

            int et = halfedges[e].twin;
            int et_n = halfedges[et].next;
            int et_p = halfedges[et_n].next;

            int f1 = halfedges[e].face;
            int f2 = halfedges[et].face;

            int A = halfedges[e_p].vert;
            int B = halfedges[e_n].vert;
            int D = halfedges[e].vert;
            int C = halfedges[et_p].vert;

            // 新增顶点 M
            Vec3 posM = (verts[B].pos + verts[D].pos) * 0.5;
            int M = (int)verts.size();
            verts.emplace_back(posM);

            // 新增两个面
            int f3 = (int)faces.size();
            faces.emplace_back();
            int f4 = (int)faces.size();
            faces.emplace_back();

            // 新增六条边
            int base = (int)halfedges.size();
            halfedges.resize(base + 6);
            int h0 = base + 0;   // D→M
            int h1 = base + 1;   // M→A
            int h2 = base + 2;   // A→M
            int h3 = base + 3;   // M→B
            int h4 = base + 4;   // M→C
            int h5 = base + 5;   // C→M

            // 修改复用的半边
            // e  原来 D→B，现在变成 B→M
            halfedges[e].vert = B;
            // et 原来 B→D，现在变成 M→D
            halfedges[et].vert = M;

            // 设置新半边的 vert
            halfedges[h0].vert = D;
            halfedges[h1].vert = M;
            halfedges[h2].vert = A;
            halfedges[h3].vert = M;
            halfedges[h4].vert = M;
            halfedges[h5].vert = C;

            // 设置 next
            // △ADM（f1）：e_p(A→D) → h0(D→M) → h1(M→A) → e_p
            halfedges[e_p].next = h0;
            halfedges[h0].next = h1;
            halfedges[h1].next = e_p;

            // △AMB（f3）：h2(A→M) → h3(M→B) → e_n(B→A) → h2
            halfedges[h2].next = h3;
            halfedges[h3].next = e_n;
            halfedges[e_n].next = h2;

            // △BMC（f2）：e(B→M) → h4(M→C) → et_p(C→B) → e
            halfedges[e].next = h4;
            halfedges[h4].next = et_p;
            halfedges[et_p].next = e;

            // △MDC（f4）：et(M→D) → et_n(D→C) → h5(C→M) → et
            halfedges[et].next = et_n;
            halfedges[et_n].next = h5;
            halfedges[h5].next = et;

            // 设置 face
            halfedges[e_p].face = f1;
            halfedges[h0].face = f1;
            halfedges[h1].face = f1;

            halfedges[h2].face = f3;
            halfedges[h3].face = f3;
            halfedges[e_n].face = f3;

            halfedges[e].face = f2;
            halfedges[h4].face = f2;
            halfedges[et_p].face = f2;

            halfedges[et].face = f4;
            halfedges[et_n].face = f4;
            halfedges[h5].face = f4;

            // 设置 twin
            halfedges[h0].twin = et;   // D→M ↔ M→D(et)
            halfedges[et].twin = h0;
            halfedges[h1].twin = h2;       // M→A ↔ A→M
            halfedges[h2].twin = h1;
            halfedges[h3].twin = e;        // M→B ↔ B→M(e)
            halfedges[e].twin = h3;
            halfedges[h4].twin = h5;       // M→C ↔ C→M
            halfedges[h5].twin = h4;

            // 更新 Face 代表半边
            faces[f1].he = e_p;
            faces[f2].he = e;
            faces[f3].he = h2;
            faces[f4].he = et;

            // 更新顶点 M 的出边
            verts[M].outHE = h3;

            // 更新顶点 B 和 D 的出边（原来可能指向 e 或 et）
            // D 的出边若原来指向 e（现在 e 起点变成 B），需要更新
            if (verts[D].outHE == e)
                verts[D].outHE = h0;    // h0 = D→M，从 D 出发

            // B 的出边若原来指向 et（现在 et 起点变成 M），需要更新
            if (verts[B].outHE == et)
                verts[B].outHE = e_n;   // e_n = B→A，从 B 出发

            return M;
        }

        // 顶点 v 的所有出边
        std::vector<int> outEdge(int v) const {
            std::vector<int> result;
            int start = verts[v].outHE;
            if (start == -1) return result;

            int rightmost = findRightmostHE(v);

            if (rightmost == -1) {
                // 内部顶点：逆时针 (CCW) 遍历
                int cur = start;
                do
                {
                    result.push_back(cur);
                    int prev = halfedges[halfedges[cur].next].next;
                    cur = halfedges[prev].twin;
                } while (cur != start);
            }
            else
            {
                // 边界顶点：从最右侧边界出边逆时针 (CCW) 遍历
                int cur = rightmost;
                while (true)
                {
                    result.push_back(cur);
                    int prev = halfedges[halfedges[cur].next].next;
                    int tw = halfedges[prev].twin;
                    if (tw == -1) break;
                    cur = tw;
                }
            }
            return result;
        }

        // 边折叠
        // 将半边 he 的起点 u 合并到终点 v，并标记相关的两个面和边为失效 (-1)
        // 返回：合并后的顶点索引 (v)；如果因为流形条件限制无法折叠，则返回 -1
        int collapseEdge(int he) {
            int e = he;
            int u = halfedges[e].vert;
            int e_n = halfedges[e].next;
            int v = halfedges[e_n].vert;

            // Step 1: Link Condition 检查，防止生成非流形网格
            std::vector<int> Nu = neighborVerts(u);
            std::vector<int> Nv = neighborVerts(v);
            int shareCount = 0;
            for (int nu : Nu) {
                for (int nv : Nv) {
                    if (nu == nv) shareCount++;
                }
            }
            // 如果是内部边，期望恰好有 2 个公共点；如果是边界边，期望 1 个
            int expectedShared = (halfedges[e].twin != -1) ? 2 : 1;
            if (shareCount != expectedShared) return -1; // 违反 Link Condition，强行折叠会导致非流形，拒绝操作

            // Step 2: 提前收集拓扑上下文 (在修改任何指针之前)
            int e_p = halfedges[e_n].next;
            int f_top = halfedges[e].face;
            int w1 = halfedges[e_p].vert;
            int tw_n = halfedges[e_n].twin;
            int tw_p = halfedges[e_p].twin;

            int et = halfedges[e].twin;
            int f_bot = -1, et_n = -1, et_p = -1, w2 = -1, tw_tn = -1, tw_tp = -1;
            if (et != -1) {
                f_bot = halfedges[et].face;
                et_n = halfedges[et].next;
                et_p = halfedges[et_n].next;
                w2 = halfedges[et_p].vert;
                tw_tn = halfedges[et_n].twin;
                tw_tp = halfedges[et_p].twin;
            }

            // 记录所有即将被删除的半边
            std::vector<int> deleted_hes = { e, e_n, e_p };
            if (et != -1) {
                deleted_hes.push_back(et);
                deleted_hes.push_back(et_n);
                deleted_hes.push_back(et_p);
            }

            // 提取受影响顶点的原始出边集合 (用于后续筛选安全出边)
            std::vector<int> u_out_edges = outEdges(u);
            std::vector<int> v_out_edges = outEdges(v);
            std::vector<int> w1_out_edges = outEdges(w1);
            std::vector<int> w2_out_edges;
            if (w2 != -1) w2_out_edges = outEdges(w2);

            // 辅助 Lambda：从出边列表中挑选一条不会被删除的边
            auto getValidOutHE = [&](const std::vector<int>& out_list) -> int {
                for (int h : out_list) {
                    if (std::find(deleted_hes.begin(), deleted_hes.end(), h) == deleted_hes.end()) {
                        return h;
                    }
                }
                return -1;
            };

            // Step 3: 更新几何位置 (采用中点策略)
            verts[v].pos = (verts[u].pos + verts[v].pos) * 0.5;

            // Step 4: 转移顶点 u 的扇面
            // 所有以 u 为起点的幸存半边，现在起点改为 v
            for (int h : u_out_edges) {
                if (std::find(deleted_hes.begin(), deleted_hes.end(), h) == deleted_hes.end()) {
                    halfedges[h].vert = v;
                }
            }

            // Step 5: 缝合撕裂的拓扑 (连接外部的 Twin)
            // 缝合顶面 (f_top) 留下的缺口
            if (tw_n != -1) halfedges[tw_n].twin = tw_p;
            if (tw_p != -1) halfedges[tw_p].twin = tw_n;

            // 缝合底面 (f_bot) 留下的缺口
            if (et != -1) {
                if (tw_tn != -1) halfedges[tw_tn].twin = tw_tp;
                if (tw_tp != -1) halfedges[tw_tp].twin = tw_tn;
            }

            // Step 6: 修复受影响顶点的 outHE 指针
            // 优先从 v 原本的出边中找；如果没有，说明 v 被完全包围，从 u 转移过来的边中找
            verts[v].outHE = getValidOutHE(v_out_edges);
            if (verts[v].outHE == -1) {
                verts[v].outHE = getValidOutHE(u_out_edges);
            }

            verts[w1].outHE = getValidOutHE(w1_out_edges);
            if (w2 != -1) verts[w2].outHE = getValidOutHE(w2_out_edges);

            // Step 7: 废弃元素标记 (Lazy Deletion)
            verts[u].outHE = -1; // 顶点 u 被逻辑删除
            faces[f_top].he = -1;
            if (f_bot != -1) faces[f_bot].he = -1;

            for (int h : deleted_hes) {
                halfedges[h].face = -1;
                halfedges[h].next = -1;
                halfedges[h].twin = -1;
                halfedges[h].vert = -1;
            }

            return v;
        }

        // 找到并返回所有边界半边的索引
        std::vector<int> findBoundaryEdges() const {
            std::vector<int> result;

            for (int i = 0; i < (int)halfedges.size(); ++i) {
                if (halfedges[i].twin == -1) {
                    result.push_back(i);
                }
            }

            return result;
        }

        // 提取所有边界环 (Boundary Loops)
        // 返回：包含多个环的二维数组，每个环是一个按连通顺序排列的顶点索引集合
        std::vector<std::vector<int>> findBoundaryLoops() const {
            std::vector<std::vector<int>> loops;
            int n_he = (int)halfedges.size();

            std::vector<bool> visited(n_he, false);

            for (int i = 0; i < n_he; ++i) {
                // 寻找一个合格的起点：
                // 1. 必须是边界半边 (twin == -1)
                // 2. 还没有被归入任何已知的环 (!visited[i])
                // 3. 确保不是被惰性删除的废弃半边 (vert != -1)
                if (halfedges[i].twin == -1 && !visited[i] && halfedges[i].vert != -1) {
                    std::vector<int> current_loop;
                    int cur = i;
                    do
                    {
                        visited[cur] = true;
                        current_loop.push_back(halfedges[cur].vert);

                        // 【核心逻辑】：寻找从当前半边终点出发的“下一条边界半边”
                        int he_n = halfedges[cur].next;
                        while (halfedges[he_n].twin != -1)
                        {
                            he_n = halfedges[halfedges[he_n].twin].next;
                        }

                        cur = he_n;
                    } while (cur != i);

                    loops.push_back(current_loop);
                }
            }

            return loops;
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
    //
    // 注意：verts 是 2D 顶点数组，mesh 内部的 Vertex.pos 是 3D 的，
    //       这里用外部 verts 做 2D 几何判断以保持精度
    inline int locatePoint_walk(
        const std::vector<Vec2>& verts,
        const HalfEdgeMesh& mesh,
        const Vec2& P,
        int startFace = 0)
    {
        int nFaces = (int)mesh.faces.size();
        int f = startFace;
        int maxIter = nFaces + 10;   // 防止死循环

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

            real x_min = points[0].x; real x_max = points[0].x; real y_min = points[0].y; real y_max = points[0].y;

            for (const auto& p : points) {
                if (p.x < x_min) x_min = p.x;
                if (p.x > x_max) x_max = p.x;
                if (p.y < y_min) y_min = p.y;
                if (p.y > y_max) y_max = p.y;
            }

            real x_mid = (x_min + x_max) * 0.5;
            real y_mid = (y_min + y_max) * 0.5;

            real dx = x_max - x_min;
            real dy = y_max - y_min;
            real delta = (dx + dy + 1.0) * 10.0;

            Vec2 v0 = { x_mid - 2.0 * delta, y_mid - delta };
            Vec2 v1 = { x_mid + 2.0 * delta, y_mid - delta };
            Vec2 v2 = { x_mid,                y_mid + 2.0 * delta };

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
        // 返回：每条边界边用 {顶点索引A, 顶点索引B} 表示，保留原始方向（逆时针）
        //
        // 实现方式：
        //   用归一化的 (min,max) 键来统计出现次数（判断是否边界）
        //   同时用有向的 (va,vb) 键来记录原始方向（用于构建新三角形）
        // -----------------------------------------------------------------------
        std::vector<std::array<int, 2>> findBoundaryEdges(const std::vector<int>& cavity) {

            // 统计 cavity 内每条边出现的次数
            // 出现1次 = 边界边；出现2次 = 内部边（两侧都在 cavity 内，需删除）
            std::map<std::array<int, 2>, int> edgeCount;
            // 保留有向边，用于提取原始方向
            std::map<std::array<int, 2>, std::array<int, 2>> edgeDirected;

            for (int fi : cavity) {
                for (int i = 0; i < 3; ++i) {
                    int va = tris[fi][i];
                    int vb = tris[fi][(i + 1) % 3];
                    std::array<int, 2> key = { std::min(va, vb), std::max(va, vb) };
                    edgeCount[key]++;
                    edgeDirected[key] = { va, vb };  // 后写入的覆盖前面的，但边界边只出现1次所以没问题
                }
            }

            std::vector<std::array<int, 2>> boundary;
            for (auto& [key, cnt] : edgeCount) {
                if (cnt == 1) boundary.push_back(edgeDirected[key]);
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

            // Step 2：找到空腔的边界边（保留了原始逆时针方向）
            std::vector<std::array<int, 2>> boundary = findBoundaryEdges(cavity);

            // Step 3：删除空腔中的三角形
            // 使用 swap-with-last 技巧，O(1) 删除每个三角形，避免 vector::erase 的 O(n) 搬移
            std::sort(cavity.begin(), cavity.end(), std::greater<int>());
            for (int fi : cavity) {
                tris[fi] = tris.back();
                tris.pop_back();
            }

            // Step 4：加入新顶点
            int newV = (int)verts.size();
            verts.push_back(P);

            // Step 5：用每条边界边和新点构成新三角形
            // 边界边的方向是 cavity 三角形的逆时针方向，
            // 从 cavity 外侧看，边界边方向需要反转，因此新三角形是 vb→va→newV
            for (auto& edge : boundary) {
                int va = edge[0], vb = edge[1];
                // 从 cavity 外部看，边界边的外侧朝向是 vb→va
                // 新三角形 (vb, va, newV) 应该是逆时针的
                // 用叉积验证，如果不是则交换
                real c = Geo::cross2D(verts[vb], verts[va], verts[newV]);
                if (c > EPS) {
                    tris.push_back({ vb, va, newV });
                }
                else {
                    tris.push_back({ va, vb, newV });
                }
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
        }
        fv.close();

        std::ofstream ft(prefix + "_tris.csv");
        ft << "i0,i1,i2\n";
        for (const auto& t : tris)
            ft << t[0] << "," << t[1] << "," << t[2] << "\n";
        ft.close();

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
        real xmin = verts[0].x, xmax = verts[0].x;
        real ymin = verts[0].y, ymax = verts[0].y;
        for (const auto& v : verts) {
            if (v.x < xmin) xmin = v.x;
            if (v.x > xmax) xmax = v.x;
            if (v.y < ymin) ymin = v.y;
            if (v.y > ymax) ymax = v.y;
        }
        const real margin = svgSize * 0.08;
        const real drawSize = svgSize - 2.0 * margin;
        real range = std::max(xmax - xmin, ymax - ymin);
        if (range < 1e-9) range = 1.0;
        // 世界坐标 → SVG 坐标（y 轴翻转）
        auto sx = [&](real x) { return margin + (x - xmin) / range * drawSize; };
        auto sy = [&](real y) { return margin + (ymax - y) / range * drawSize; };

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
            real ax = sx(verts[t[0]].x), ay = sy(verts[t[0]].y);
            real bx = sx(verts[t[1]].x), by = sy(verts[t[1]].y);
            real cx = sx(verts[t[2]].x), cy = sy(verts[t[2]].y);
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