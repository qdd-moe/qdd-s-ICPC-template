## Geometry

### 2D Geometry Basics

```cpp
#define y1 qwq

using ld = double;

const ld PI = acos(-1);
const ld EPS = 1e-8;

int sgn(ld x) { return x < -EPS ? -1 : x > EPS; }

// don't use sgn directly
bool eq(ld x, ld y) { return sgn(x - y) == 0; }
bool ne(ld x, ld y) { return sgn(x - y) != 0; }
bool lt(ld x, ld y) { return sgn(x - y) < 0; }
bool gt(ld x, ld y) { return sgn(x - y) > 0; }
bool le(ld x, ld y) { return sgn(x - y) <= 0; }
bool ge(ld x, ld y) { return sgn(x - y) >= 0; }

struct V {
  ld x, y;
  constexpr V(ld x = 0, ld y = 0) : x(x), y(y) {}
  V operator+(V b) const { return V(x + b.x, y + b.y); }
  V operator-(V b) const { return V(x - b.x, y - b.y); }
  V operator*(ld k) const { return V(x * k, y * k); }
  V operator/(ld k) const { return V(x / k, y / k); }
  ld len() const { return hypot(x, y); }
  ld len2() const { return x * x + y * y; }
};

ostream& operator<<(ostream& os, V p) { return os << "(" << p.x << ", " << p.y << ")"; }
istream& operator>>(istream& is, V& p) { return is >> p.x >> p.y; }

ld dist(V a, V b) { return (b - a).len(); }
ld dot(V a, V b) { return a.x * b.x + a.y * b.y; }
ld det(V a, V b) { return a.x * b.y - a.y * b.x; }
ld cross(V s, V t, V o) { return det(s - o, t - o); }

ld to_rad(ld deg) { return deg / 180 * PI; }

// quadrant
int quad(V p) {
  int x = sgn(p.x), y = sgn(p.y);
  if (x > 0 && y >= 0) return 1;
  if (x <= 0 && y > 0) return 2;
  if (x < 0 && y <= 0) return 3;
  if (x >= 0 && y < 0) return 4;
  assert(0);
}

// sorting by polar angle
struct cmp_angle {
  V p;
  cmp_angle(V p = V()) : p(p) {}
  bool operator () (V a, V b) const {
    int qa = quad(a - p), qb = quad(b - p);
    if (qa != qb) return qa < qb;
    int d = sgn(cross(a, b, p));
    if (d) return d > 0;
    return dist(a, p) < dist(b, p);
  }
};

// unit vector
V unit(V p) { return eq(p.len(), 0) ? V(1, 0) : p / p.len(); }

// rotate conterclockwise by r radians
V rot(V p, ld r) {
  return V(p.x * cos(r) - p.y * sin(r), p.x * sin(r) + p.y * cos(r));
}
V rot_ccw90(V p) { return V(-p.y, p.x); }
V rot_cw90(V p) { return V(p.y, -p.x); }

// point on segment, le(dot(...) , 0) contains endpoints, lt(dot(...) , 0) otherwise
bool p_on_seg(V p, V a, V b) {
  return eq(det(p - a, b - a), 0) && le(dot(p - a, p - b), 0);
}

// point on ray, ge(dot(...) , 0) contains endpoints, gt(dot(...) , 0) otherwise
bool p_on_ray(V p, V a, V b) {
  return eq(det(p - a, b - a), 0) && ge(dot(p - a, b - a), 0);
}

// intersection of lines
V intersect(V a, V b, V c, V d) {
  ld s1 = cross(c, d, a), s2 = cross(c, d, b);
  return (a * s2 - b * s1) / (s2 - s1);
}

// projection of a point onto a line
V proj(V p, V a, V b) {
  return a + (b - a) * dot(b - a, p - a) / (b - a).len2();
}

// symmetric point about a line
V reflect(V p, V a, V b) {
  return proj(p, a, b) * 2 - p;
}

// closest point on segment
V closest_point_on_seg(V p, V a, V b) {
  if (lt(dot(b - a, p - a), 0)) return a;
  if (lt(dot(a - b, p - b), 0)) return b;
  return proj(p, a, b);
}

// centroid
V centroid(V a, V b, V c) {
  return (a + b + c) / 3;
}

// incenter
V incenter(V a, V b, V c) {
  ld AB = dist(a, b), AC = dist(a, c), BC = dist(b, c);
  // ld r = abs(cross(b, c, a)) / (AB + AC + BC);
  return (a * BC + b * AC + c * AB) / (AB + BC + AC);
}

// circumcenter
V circumcenter(V a, V b, V c) {
  V mid1 = (a + b) / 2, mid2 = (a + c) / 2;
  // ld r = dist(a, b) * dist(b, c) * dist(c, a) / 2 / abs(cross(b, c, a));
  return intersect(mid1, mid1 + rot_ccw90(b - a), mid2, mid2 + rot_ccw90(c - a));
}

// orthocenter
V orthocenter(V a, V b, V c) {
  return centroid(a, b, c) * 3 - circumcenter(a, b, c) * 2;
}

// excenters (opposite to a, b, c)
vector<V> excenters(V a, V b, V c) {
  ld AB = dist(a, b), AC = dist(a, c), BC = dist(b, c);
  V p1 = (a * (-BC) + b * AC + c * AB) / (AB + AC - BC);
  V p2 = (a * BC + b * (-AC) + c * AB) / (AB - AC + BC);
  V p3 = (a * BC + b * AC + c * (-AB)) / (-AB + AC + BC);
  return {p1, p2, p3};
}
```

### Polygons

```cpp
// polygon area
ld area(const vector<V>& s) {
  ld ret = 0;
  for (int i = 0; i < s.size(); i++) {
    ret += det(s[i], s[(i + 1) % s.size()]);
  }
  return ret / 2;
}

// polygon centroid
V centroid(const vector<V>& s) {
  V c;
  for (int i = 0; i < s.size(); i++) {
    c = c + (s[i] + s[(i + 1) % s.size()]) * det(s[i], s[(i + 1) % s.size()]);
  }
  return c / 6.0 / area(s);
}

// point and polygon
// 1 inside 0 on border -1 outside
int inside(const vector<V>& s, V p) {
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

// convex hull, points cannot be duplicated
// lt(cross(...), 0) allow point on edges le(cross(...), 0) otherwise
// will change the order of the input points
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

// is convex?
bool is_convex(const vector<V>& s) {
  for (int i = 0; i < s.size(); i++) {
    if (lt(cross(s[(i + 1) % s.size()], s[(i + 2) % s.size()], s[i]), 0)) return false;
  }
  return true;
}

// point and convex hull
// 1 inside 0 on border -1 outside
int inside(const vector<V>& s, V p) {
  for (int i = 0; i < s.size(); i++) {
    if (lt(cross(s[i], s[(i + 1) % s.size()], p), 0)) return -1;
    if (p_on_seg(p, s[i], s[(i + 1) % s.size()])) return 0;
  }
  return 1;
}

// closest pair of points, sort by x-coordinate first
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
  sort(q.begin(), q.end(), [](auto& a, auto& b) { return a.y < b.y; });
  for (int i = 1; i < q.size(); i++) {
    for (int j = i - 1; j >= 0 && q[j].y >= q[i].y - ret; j--) {
      ret = min(ret, dist(q[i], q[j]));
    }
  }
  return ret;
}
```

### Circles

```cpp
struct C {
  V o;
  ld r;
  C(V o, ld r) : o(o), r(r) {}
};

// sector area, radius r, angle d
ld area_sector(ld r, ld d) { return r * r * d / 2; }

// find the tangent line to a circle and return the tangent point
vector<V> tangent_point(C c, V p) {
  ld k = c.r / dist(c.o, p);
  if (gt(k, 1)) return vector<V>();
  if (eq(k, 1)) return {p};
  V a = (p - c.o) * k;
  return {c.o + rot(a, acos(k)), c.o + rot(a, -acos(k))};
}

// minimum covering circle
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

### 3D Geometry

```cpp
struct V {
  ld x, y, z;
  constexpr V(ld x = 0, ld y = 0, ld z = 0) : x(x), y(y), z(z) {}
  V operator+(V b) const { return V(x + b.x, y + b.y, z + b.z); }
  V operator-(V b) const { return V(x - b.x, y - b.y, z - b.z); }
  V operator*(ld k) const { return V(x * k, y * k, z * k); }
  V operator/(ld k) const { return V(x / k, y / k, z / k); }
  ld len() const { return sqrt(len2()); }
  ld len2() const { return x * x + y * y + z * z; }
};

ostream& operator<<(ostream& os, V p) { return os << "(" << p.x << "," << p.y << "," << p.z << ")"; }
istream& operator>>(istream& is, V& p) { return is >> p.x >> p.y >> p.z; }

ld dist(V a, V b) { return (b - a).len(); }
ld dot(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
V det(V a, V b) { return V(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }
V cross(V s, V t, V o) { return det(s - o, t - o); }
ld mix(V a, V b, V c) { return dot(a, det(b, c)); }
```
