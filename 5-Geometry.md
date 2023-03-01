## 计算几何

### 二维几何基础

```cpp
#define y1 qwq

using ld = double;

const ld PI = acos(-1);
const ld EPS = 1e-8;

int sgn(ld x) { return x < -EPS ? -1 : x > EPS; }

// 不要直接使用sgn
bool eq(ld x, ld y) { return sgn(x - y) == 0; }
bool ne(ld x, ld y) { return sgn(x - y) != 0; }
bool lt(ld x, ld y) { return sgn(x - y) < 0; }
bool gt(ld x, ld y) { return sgn(x - y) > 0; }
bool le(ld x, ld y) { return sgn(x - y) <= 0; }
bool ge(ld x, ld y) { return sgn(x - y) >= 0; }

struct V {
  ld x, y;
  constexpr V(ld x = 0, ld y = 0) : x(x), y(y) {}
  V operator+(const V& b) const { return V(x + b.x, y + b.y); }
  V operator-(const V& b) const { return V(x - b.x, y - b.y); }
  V operator*(ld k) const { return V(x * k, y * k); }
  V operator/(ld k) const { return V(x / k, y / k); }
  ld len() const { return hypot(x, y); }
  ld len2() const { return x * x + y * y; }
};

ostream& operator<<(ostream& os, const V& p) { return os << "(" << p.x << ", " << p.y << ")"; }
istream& operator>>(istream& is, V& p) { return is >> p.x >> p.y; }

ld dist(const V& a, const V& b) { return (b - a).len(); }
ld dot(const V& a, const V& b) { return a.x * b.x + a.y * b.y; }
ld det(const V& a, const V& b) { return a.x * b.y - a.y * b.x; }
ld cross(const V& s, const V& t, const V& o) { return det(s - o, t - o); }

ld to_rad(ld deg) { return deg / 180 * PI; }

// 象限
int quad(const V& p) {
  int x = sgn(p.x), y = sgn(p.y);
  if (x > 0 && y >= 0) return 1;
  if (x <= 0 && y > 0) return 2;
  if (x < 0 && y <= 0) return 3;
  if (x >= 0 && y < 0) return 4;
  assert(0);
}

// 极角排序
struct cmp_angle {
  V p;
  cmp_angle(const V& p = V()) : p(p) {}
  bool operator () (const V& a, const V& b) const {
    int qa = quad(a - p), qb = quad(b - p);
    if (qa != qb) return qa < qb;
    int d = sgn(cross(a, b, p));
    if (d) return d > 0;
    return dist(a, p) < dist(b, p);
  }
};

// 单位向量
V unit(const V& p) { return eq(p.len(), 0) ? V(1, 0) : p / p.len(); }

// 逆时针旋转 r 弧度
V rot(const V& p, ld r) {
  return V(p.x * cos(r) - p.y * sin(r), p.x * sin(r) + p.y * cos(r));
}
V rot_ccw90(const V& p) { return V(-p.y, p.x); }
V rot_cw90(const V& p) { return V(p.y, -p.x); }

// 点在线段上 le(dot(...), 0) 包含端点 lt(dot(...), 0) 则不包含
bool p_on_seg(const V& p, const V& a, const V& b) {
  return eq(det(p - a, b - a), 0) && le(dot(p - a, p - b), 0);
}

// 点在射线上 ge(dot(...), 0) 包含端点 gt(dot(...), 0) 则不包含
bool p_on_ray(const V& p, const V& a, const V& b) {
  return eq(det(p - a, b - a), 0) && ge(dot(p - a, b - a), 0);
}

// 求直线交点
V intersect(const V& a, const V& b, const V& c, const V& d) {
  ld s1 = cross(c, d, a), s2 = cross(c, d, b);
  return (a * s2 - b * s1) / (s2 - s1);
}

// 点在直线上的投影点
V proj(const V& p, const V& a, const V& b) {
  return a + (b - a) * dot(b - a, p - a) / (b - a).len2();
}

// 点关于直线的对称点
V reflect(const V& p, const V& a, const V& b) {
  return proj(p, a, b) * 2 - p;
}

// 点到线段的最近点
V closest_point_on_seg(const V& p, const V& a, const V& b) {
  if (lt(dot(b - a, p - a), 0)) return a;
  if (lt(dot(a - b, p - b), 0)) return b;
  return proj(p, a, b);
}

// 三角形重心
V centroid(const V& a, const V& b, const V& c) {
  return (a + b + c) / 3;
}

// 内心
V incenter(const V& a, const V& b, const V& c) {
  ld AB = dist(a, b), AC = dist(a, c), BC = dist(b, c);
  // ld r = abs(cross(b, c, a)) / (AB + AC + BC);
  return (a * BC + b * AC + c * AB) / (AB + BC + AC);
}

// 外心
V circumcenter(const V& a, const V& b, const V& c) {
  V mid1 = (a + b) / 2, mid2 = (a + c) / 2;
  // ld r = dist(a, b) * dist(b, c) * dist(c, a) / 2 / abs(cross(b, c, a));
  return intersect(mid1, mid1 + rot_ccw90(b - a), mid2, mid2 + rot_ccw90(c - a));
}

// 垂心
V orthocenter(const V& a, const V& b, const V& c) {
  return centroid(a, b, c) * 3 - circumcenter(a, b, c) * 2;
}

// 旁心（三个）
vector<V> escenter(const V& a, const V& b, const V& c) {
  ld AB = dist(a, b), AC = dist(a, c), BC = dist(b, c);
  V p1 = (a * (-BC) + b * AC + c * AB) / (AB + AC - BC);
  V p2 = (a * BC + b * (-AC) + c * AB) / (AB - AC + BC);
  V p3 = (a * BC + b * AC + c * (-AB)) / (-AB + AC + BC);
  return {p1, p2, p3};
}
```

### 多边形

```cpp
// 多边形面积
ld area(const vector<V>& s) {
  ld ret = 0;
  for (int i = 0; i < s.size(); i++) {
    ret += det(s[i], s[(i + 1) % s.size()]);
  }
  return ret / 2;
}

// 多边形重心
V centroid(const vector<V>& s) {
  V c;
  for (int i = 0; i < s.size(); i++) {
    c = c + (s[i] + s[(i + 1) % s.size()]) * det(s[i], s[(i + 1) % s.size()]);
  }
  return c / 6.0 / area(s);
}

// 点是否在多边形中
// 1 inside 0 on border -1 outside
int inside(const vector<V>& s, const V& p) {
  int cnt = 0;
  for (int i = 0; i < s.size(); i++) {
    V a = s[i], b = s[(i + 1) % s.size()];
    if (p_on_seg(p, a, b)) return 0;
    if (le(a.y, b.y)) swap(a, b);
    if (gt(p.y, a.y)) continue;
    if (le(p.y, b.y)) continue;
    cnt += gt(cross(b, a, p), 0);
  }
  return (cnt & 1) ? 1 : -1;
}

// 构建凸包 点不可以重复
// lt(cross(...), 0) 边上可以有点 le(cross(...), 0) 则不能
// 会改变输入点的顺序
vector<V> convex_hull(vector<V>& s) {
  // assert(s.size() >= 3);
  sort(s.begin(), s.end(), [](V &a, V &b) { return eq(a.x, b.x) ? lt(a.y, b.y) : lt(a.x, b.x); });
  vector<V> ret(2 * s.size());
  int sz = 0;
  for (int i = 0; i < s.size(); i++) {
    while (sz > 1 && le(cross(ret[sz - 1], s[i], ret[sz - 2]), 0)) sz--;
    ret[sz++] = s[i];
  }
  int k = sz;
  for (int i = s.size() - 2; i >= 0; i--) {
    while (sz > k && le(cross(ret[sz - 1], s[i], ret[sz - 2]), 0)) sz--;
    ret[sz++] = s[i];
  }
  ret.resize(sz - (s.size() > 1));
  return ret;
}

// 多边形是否为凸包
bool is_convex(const vector<V>& s) {
  for (int i = 0; i < s.size(); i++) {
    if (lt(cross(s[(i + 1) % s.size()], s[(i + 2) % s.size()], s[i]), 0)) return false;
  }
  return true;
}

// 点是否在凸包中
// 1 inside 0 on border -1 outside
int inside(const vector<V>& s, const V& p) {
  for (int i = 0; i < s.size(); i++) {
    if (lt(cross(s[i], s[(i + 1) % s.size()], p), 0)) return -1;
    if (p_on_seg(p, s[i], s[(i + 1) % s.size()])) return 0;
  }
  return 1;
}

// 最近点对, 先要按照 x 坐标排序
// min_dist(s, 0, s.size())
ld min_dist(const vector<V>& s, int l, int r) {
  if (r - l <= 5) {
    ld ret = 1e100;
    for (int i = l; i < r; i++) {
      for (int j = i + 1; j < r; j++) {
        ret = min(ret, dist(s[i], s[j]));
      }
    }
    return ret;
  }
  int m = (l + r) >> 1;
  ld ret = min(min_dist(s, l, m), min_dist(s, m, r));
  vector<V> q;
  for (int i = l; i < r; i++) {
    if (abs(s[i].x - s[m].x) <= ret) q.push_back(s[i]);
  }
  sort(q.begin(), q.end(), [](const V& a, const V& b) { return a.y < b.y; });
  for (int i = 1; i < q.size(); i++) {
    for (int j = i - 1; j >= 0 && q[j].y >= q[i].y - ret; j--) {
      ret = min(ret, dist(q[i], q[j]));
    }
  }
  return ret;
}
```

### 圆

```cpp
struct C {
  V o;
  ld r;
  C(const V& o, ld r) : o(o), r(r) {}
};

// 扇形面积，半径 r 圆心角 d
ld area_sector(ld r, ld d) { return r * r * d / 2; }

// 过一点求圆的切线，返回切点
vector<V> tangent_point(const C& c, const V& p) {
  ld k = c.r / dist(c.o, p);
  if (gt(k, 1)) return vector<V>();
  if (eq(k, 1)) return {p};
  V a = (p - c.o) * k;
  return {c.o + rot(a, acos(k)), c.o + rot(a, -acos(k))};
}

// 最小圆覆盖
C min_circle_cover(vector<V> a) {
  shuffle(a.begin(), a.end(), rng);
  V o = a[0];
  ld r = 0;
  int n = a.size();
  for (int i = 1; i < n; i++) if (gt(dist(a[i], o), r)) {
    o = a[i]; r = 0;
    for (int j = 0; j < i; j++) if (gt(dist(a[j], o), r)) {
      o = (a[i] + a[j]) / 2;
      r = dist(a[j], o);
      for (int k = 0; k < j; k++) if (gt(dist(a[k], o), r)) {
        o = circumcenter(a[i], a[j], a[k]);
        r = dist(a[k], o);
      }
    }
  }
  return C(o, r);
}
```

### 三维几何

```cpp
struct V {
  ld x, y, z;
  constexpr V(ld x = 0, ld y = 0, ld z = 0) : x(x), y(y), z(z) {}
  V operator+(const V& b) const { return V(x + b.x, y + b.y, z + b.z); }
  V operator-(const V& b) const { return V(x - b.x, y - b.y, z - b.z); }
  V operator*(ld k) const { return V(x * k, y * k, z * k); }
  V operator/(ld k) const { return V(x / k, y / k, z / k); }
  ld len() const { return sqrt(len2()); }
  ld len2() const { return x * x + y * y + z * z; }
};

ostream& operator<<(ostream& os, const V& p) { return os << "(" << p.x << "," << p.y << "," << p.z << ")"; }
istream& operator>>(istream& is, V& p) { return is >> p.x >> p.y >> p.z; }

ld dist(const V& a, const V& b) { return (b - a).len(); }
ld dot(const V& a, const V& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
V det(const V& a, const V& b) { return V(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }
V cross(const V& s, const V& t, const V& o) { return det(s - o, t - o); }
ld mix(const V& a, const V& b, const V& c) { return dot(a, det(b, c)); }
```
