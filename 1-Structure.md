## 数据结构

**Unless specified, all data structure interfaces are 0-indexed, [l, r]**

### 并查集

```cpp
struct dsu {
  vector<int> p;
  dsu(int n) : p(n, -1) {}
  int find(int x) { return (p[x] < 0) ? x : p[x] = find(p[x]); }
  bool merge(int x, int y) {
    x = find(x), y = find(y);
    if (x == y) return 0;
    p[y] += p[x];
    p[x] = y;
    return 1;
  }
};
```

+ 动态开点并查集

```cpp
unordered_map<int, int> pa;

void _set(int x) { if (!pa.count(x)) pa[x] = -1; }
int find(int x) { return (pa[x] < 0) ? x : pa[x] = find(pa[x]); }

void merge(int x, int y) {
  x = find(x), y = find(y);
  if (x == y) return;
  if (pa[x] < pa[y]) swap(x, y);
  pa[y] += pa[x];
  pa[x] = y;
}
```

+ 可回滚并查集

```cpp
struct dsu {
  vector<int> p, sz, undo;
  dsu(int n) : p(n, -1), sz(n, 1) {}
  int find(int x) { return (p[x] < 0) ? x : find(p[x]); }
  bool merge(int x, int y) {
    x = find(x), y = find(y);
    if (x == y) return 0;
    if (sz[x] > sz[y]) swap(x, y);
    undo.push_back(x);
    p[x] = y;
    sz[y] += sz[x];
    return 1;
  }
  void rollback() {
    int x = undo.back();
    undo.pop_back();
    sz[p[x]] -= sz[x];
    p[x] = -1;
  }
};
```

### RMQ

+ 一维

```cpp
struct RMQ {
  vector<vector<int>> st;

  RMQ(const vector<int>& a) {
    int n = a.size();
    st.assign(__lg(n) + 1, vector<int>(n));
    for (int i = 0; i < n; i++) {
      st[0][i] = a[i];
    }
    for (int k = 1; (1 << k) <= n; k++) {
      for (int i = 0; i + (1 << k) <= n; i++) {
        st[k][i] = max(st[k - 1][i], st[k - 1][i + (1 << (k - 1))]);
      }
    }
  }

  int query(int l, int r) {
    int k = __lg(r - l + 1);
    return max(st[k][l], st[k][r - (1 << k) + 1]);
  }
};
```

+ 二维

```cpp
struct RMQ {
  int st[11][11][N][N]; // 11 = ((int)log2(N) + 1)

  void init(int n, int m) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        st[0][0][i][j] = a[i][j];
      }
    }
    for (int i = 0; (1 << i) <= n; i++) {
      for (int j = 0; (1 << j) <= m; j++) {
        if (i == 0 && j == 0) continue;
        for (int r = 0; r + (1 << i) - 1 < n; r++) {
          for (int c = 0; c + (1 << j) - 1 < m; c++) {
            if (i == 0) {
              st[i][j][r][c] = max(st[i][j - 1][r][c], st[i][j - 1][r][c + (1 << (j - 1))]);
            } else {
              st[i][j][r][c] = max(st[i - 1][j][r][c], st[i - 1][j][r + (1 << (i - 1))][c]);
            }
          }
        }
      }
    }
  }

  int query(int r1, int c1, int r2, int c2) {
    int x = __lg(r2 - r1 + 1);
    int y = __lg(c2 - c1 + 1);
    int m1 = st[x][y][r1][c1];
    int m2 = st[x][y][r1][c2 - (1 << y) + 1];
    int m3 = st[x][y][r2 - (1 << x) + 1][c1];
    int m4 = st[x][y][r2 - (1 << x) + 1][c2 - (1 << y) + 1];
    return max({m1, m2, m3, m4});
  }
};
```

+ 滑动窗口 RMQ

```cpp
// k 为滑动窗口的大小
deque<int> q;
for (int i = 0, j = 0; i + k <= n; i++) {
  while (j < i + k) {
    while (!q.empty() && a[q.back()] < a[j]) q.pop_back(); // 最小值取'>'号
    q.push_back(j++);
  }
  while (q.front() < i) q.pop_front();
  rmq.push_back(a[q.front()]);
}
```

### 树状数组

+ 单点修改，区间和

```cpp
template <class T>
struct fenwick {
  int n;
  vector<T> t;
  fenwick(int n) : n(n), t(n) {}
  void add(int p, T x) {
    // assert(0 <= p && p < n);
    for (; p < n; p |= p + 1) t[p] += x;
  }
  T sum(int p) {
    T ans = 0;
    for (; p >= 0; p = (p & (p + 1)) - 1) ans += t[p];
    return ans;
  }
  void set(int p, T x) { add(p, x - range_sum(p, p)); }
  T range_sum(int l, int r) { return sum(r) - sum(l - 1); }

  // smallest is 0-th
  int kth(T k) {
    int p = 0;
    for (int i = 1 << __lg(n); i; i /= 2) {
      if (p + i <= n && k >= t[p + i - 1]) {
        p += i;
        k -= t[p - 1];
      }
    }
    return p;
  }
};
```

+ 区间加，单点查询

```cpp
void range_add(int l, int r, i64 x) {
  add(l, x);
  add(r + 1, -x);
}
```

+ 区间加，区间和

```cpp
fenwick t1, t2;

void range_add(int l, int r, i64 x) {
  t1.add(l, x);
  t2.add(l, l * x);
  t1.add(r + 1, -x);
  t2.add(r + 1, (r + 1) * -x);
}

i64 range_sum(int l, int r) {
  return (r + 1) * t1.sum(r) - t2.sum(r) - l * t1.sum(l - 1) + t2.sum(l - 1);
}
```

+ 二维

```cpp
template <class T>
struct fenwick_2d {
  int n, m;
  vector<vector<T>> t;

  fenwick_2d(int n, int m) : n(n), m(m), t(n, vector<T>(m)) {}

  void add(int x, int y, T d) {
    for (int i = x; i < n; i |= i + 1)
      for (int j = y; j < m; j |= j + 1) t[i][j] += d;
  }

  T sum(int x, int y) {
    T ans = 0;
    for (int i = x; i >= 0; i = (i & (i + 1)) - 1)
      for (int j = y; j >= 0; j = (j & (j + 1)) - 1) ans += t[i][j];
    return ans;
  }

  T range_sum(int x, int y, int xx, int yy) {
    return sum(xx, yy) - sum(x - 1, yy) - sum(xx, y - 1) + sum(x - 1, y - 1);
  }
};
```

+ 二维区间加，区间和

```cpp
fenwick t0, t1, t2, t3;

void add4(int x, int y, i64 d) {
  t0.add(x, y, d);
  t1.add(x, y, d * x);
  t2.add(x, y, d * y);
  t3.add(x, y, d * x * y);
}

void range_add(int x, int y, int xx, int yy, i64 d) {
  add4(x, y, d);
  add4(x, yy + 1, -d);
  add4(xx + 1, y, -d);
  add4(xx + 1, yy + 1, d);
}

i64 sum4(int x, int y) {
  return (x + 1) * (y + 1) * t0.sum(x, y)
  - (y + 1) * t1.sum(x, y)
  - (x + 1) * t2.sum(x, y)
  + t3.sum(x, y);
}

i64 range_sum(int x, int y, int xx, int yy) {
  return sum4(xx, yy) - sum4(x - 1, yy) - sum4(xx, y - 1) + sum4(x - 1, y - 1);
}
```

### 线段树

+ zkw 线段树

```cpp
struct segt {
  int n;
  vector<int> d;
  segt(int n) : n(n), d(2 * n, (int)-2e9) {}

  void set(int p, int x) {
    for (d[p += n] = x; p > 1; p /= 2) {
      d[p / 2] = max(d[p], d[p ^ 1]);
    }
  }

  int query(int l, int r) {
    int res = -2e9;
    for (l += n, r += n + 1; l < r; l /= 2, r /= 2) {
      if (l & 1) res = max(res, d[l++]);
      if (r & 1) res = max(res, d[--r]);
    }
    return res;
  }
};
```

+ 单点修改，RMQ

```cpp
// 必要时使用 const S& 卡常数
template <class S>
struct segtree {
  using OP = S (*)(S, S);
  using E = S (*)();

  int _n, size;
  vector<S> d;

  const OP op;
  const E e;

  void pull(int p) { d[p] = op(d[p * 2], d[p * 2 + 1]); }

  segtree(int n, OP op, E e) : _n(n), op(op), e(e) {
    size = 1;
    while (size < n) size <<= 1;
    d.assign(2 * size, e());
  }

  template <class I>
  segtree(I first, I last, OP op, E e) : segtree(last - first, op, e) {
    copy(first, last, d.begin() + size);
    for (int i = size - 1; i >= 1; i--) pull(i);
  }

  S ask(int ql, int qr, int p, int l, int r) {
    if (ql > r || qr < l) return e();
    if (ql <= l && qr >= r) return d[p];
    S vl = ask(ql, qr, p * 2, l, (l + r) >> 1);
    S vr = ask(ql, qr, p * 2 + 1, ((l + r) >> 1) + 1, r);
    return op(vl, vr);
  }

  void set(int p, S x) {
    p += size;
    d[p] = x;
    for (p >>= 1; p > 0; p >>= 1) pull(p);
  }

  S get(int p) { return d[p + size]; }
  S query(int l, int r) { return ask(l, r, 1, 0, size - 1); }

  // f(e()) = false
  // find the smallest r such that f(sum([l...r])) = true
  template <class F>
  int find_right(int l, F f) {
    l += size;
    S s = e();
    do {
      while (l % 2 == 0) l >>= 1;
      if (f(op(s, d[l]))) {
        while (l < size) {
          l *= 2;
          if (!f(op(s, d[l]))) {
            s = op(s, d[l]);
            l++;
          }
        }
        return l - size;
      }
      s = op(s, d[l]);
      l++;
    } while ((l & -l) != l);
    return _n;
  }

  // find the largest l such that f(sum([l...r])) = true
  template <class F>
  int find_left(int r, F f) {
    r += size + 1;
    S s = e();
    do {
      r--;
      while (r > 1 && (r % 2)) r >>= 1;
      if (f(op(d[r], s))) {
        while (r < size) {
          r = 2 * r + 1;
          if (!f(op(d[r], s))) {
            s = op(d[r], s);
            r--;
          }
        }
        return r - size;
      }
      s = op(d[r], s);
    } while ((r & -r) != r);
    return -1;
  }
};

auto op = [](int a, int b) { return max(a, b); };
auto e = []() { return -INF; };

segtree<int> seg(n, op, e);
```

+ 区间加，区间和

```cpp
template <class S, class F>
struct lazy_segtree {
#define args int p, int l, int r
#define lc p * 2, l, (l + r) >> 1
#define rc p * 2 + 1, ((l + r) >> 1) + 1, r

  using OP = S (*)(S, S);
  using E = S (*)();
  using MAP = S (*)(S, F, int, int);
  using COM = F (*)(F, F);
  using ID = F (*)();

  int _n, size, log;
  vector<S> d;
  vector<F> lz;

  const OP op;
  const E e;
  const MAP mapping;
  const COM composition;
  const ID id;

  void pull(int p) { d[p] = op(d[p * 2], d[p * 2 + 1]); }

  void apply(F f, args) {
    d[p] = mapping(d[p], f, l, r);
    if (p < size) lz[p] = composition(f, lz[p]);
  }

  void push(args) {
    if (lz[p] == id()) return;
    apply(lz[p], lc);
    apply(lz[p], rc);
    lz[p] = id();
  }

  S ask(int ql, int qr, args) {
    if (ql > r || qr < l) return e();
    if (ql <= l && qr >= r) return d[p];
    push(p, l, r);
    S vl = ask(ql, qr, lc);
    S vr = ask(ql, qr, rc);
    return op(vl, vr);
  }

  void modify(int ql, int qr, F f, args) {
    if (ql > r || qr < l) return;
    if (ql <= l && qr >= r) {
      apply(f, p, l, r);
      return;
    }
    push(p, l, r);
    modify(ql, qr, f, lc);
    modify(ql, qr, f, rc);
    pull(p);
  }

#undef args
#undef lc
#undef rc

  lazy_segtree(int n, OP op, E e, MAP mapping, COM composition, ID id)
      : _n(n), op(op), e(e), mapping(mapping), composition(composition), id(id) {
    size = 1, log = 0;
    while (size < n) size <<= 1, log++;
    d.assign(2 * size, e());
    lz.assign(size, id());
  }

  template <class I>
  lazy_segtree(I first, I last, OP op, E e, MAP mapping, COM composition, ID id)
      : lazy_segtree(last - first, op, e, mapping, composition, id) {
    copy(first, last, d.begin() + size);
    for (int i = size - 1; i >= 1; i--) pull(i);
  }

  void update(int l, int r, F f) { modify(l, r, f, 1, 0, size - 1); }
  S query(int l, int r) { return ask(l, r, 1, 0, size - 1); }

  // for binary search
  void push(int p) {
    int x = __lg(p), y = size >> x, z = p - (1 << x);
    push(p, y * z + 1, y * (z + 1));
  }

  // f(e()) = false
  // find the smallest r such that f(sum([l...r])) = true
  template <class G>
  int find_right(int l, G f) {
    l += size;
    for (int i = log; i >= 1; i--) push(l >> i);
    S s = e();
    do {
      while (l % 2 == 0) l >>= 1;
      if (f(op(s, d[l]))) {
        while (l < size) {
          push(l);
          l *= 2;
          if (!f(op(s, d[l]))) {
            s = op(s, d[l]);
            l++;
          }
        }
        return l - size;
      }
      s = op(s, d[l]);
      l++;
    } while ((l & -l) != l);
    return _n;
  }

  // find the largest l such that f(sum([l...r])) = true
  template <class G>
  int find_left(int r, G f) {
    r += size + 1;
    for (int i = log; i >= 1; i--) push((r - 1) >> i);
    S s = e();
    do {
      r--;
      while (r > 1 && (r % 2)) r >>= 1;
      if (f(op(d[r], s))) {
        while (r < size) {
          push(r);
          r = 2 * r + 1;
          if (!f(op(d[r], s))) {
            s = op(d[r], s);
            r--;
          }
        }
        return r - size;
      }
      s = op(d[r], s);
    } while ((r & -r) != r);
    return -1;
  }
};

auto op = [](i64 a, i64 b) { return a + b; };
auto e = []() { return 0LL; };
auto mapping = [](i64 a, i64 f, int l, int r) { return a + f * (r - l + 1); };
auto composition = [](i64 f, i64 g) { return f + g; };  // f：外层函数 g：内层函数
auto id = []() { return 0LL; };  // 如果是区间赋值，选取一个数据范围外的值

lazy_segtree<i64, i64> seg(n, op, e, mapping, composition, id);
```

+ 区间乘混加，区间和取模

```cpp
// S = i64, F = pair<i64, i64>
S op(S a, S b) { return (a + b) % P; }
S e() { return 0; }
S mapping(S a, F f, int l, int r) { return (a * f.first % P + (r - l + 1) * f.second % P) % P; }
F composition(F f, F g) { return F{(g.first * f.first) % P, (g.second * f.first % P + f.second) % P}; }
F id() { return F{1, 0}; }
```

### 动态开点线段树

```cpp
struct Node {
  int lc, rc, val;
  Node(int lc = 0, int rc = 0, int val = 0) : lc(lc), rc(rc), val(val) {}
} t[20 * N];

int cnt;

struct SegT {
#define mid ((pl + pr) >> 1)

  int rt, size;

  SegT(int sz) : rt(0) {
    size = 1;
    while (size < sz) size <<= 1;
  }

  int modify(int p, int pl, int pr, int k, int val) {
    if (pl > k || pr < k) return p;
    if (!p) p = ++cnt;
    if (pl == pr) t[p].val = val;
    else {
      t[p].lc = modify(t[p].lc, pl, mid, k, val);
      t[p].rc = modify(t[p].rc, mid + 1, pr, k, val);
      t[p].val = max(t[t[p].lc].val, t[t[p].rc].val);
    }
    return p;
  }

  int ask(int p, int pl, int pr, int l, int r) {
    if (l > pr || r < pl) return -INF;
    if (l <= pl && r >= pr) return t[p].val;
    int vl = ask(t[p].lc, pl, mid, l, r);
    int vr = ask(t[p].rc, mid + 1, pr, l, r);
    return max(vl, vr);
  }

  void update(int k, int val) { rt = modify(rt, 0, size - 1, k, val); }
  int query(int l, int r) { return ask(rt, 0, size - 1, l, r); }

#undef mid
};
```

### 主席树

```cpp
struct Node {
  int lc, rc, val;
  Node(int lc = 0, int rc = 0, int val = 0) : lc(lc), rc(rc), val(val) {}
} t[40 * MAXN]; // (4 + log(size)) * MAXN 小心 MLE

int cnt;

struct FST {
#define mid ((pl + pr) >> 1)

  int size;
  vector<int> root;

  FST(int sz) {
    size = 1;
    while (size < sz) size <<= 1;
    root.push_back(N(0, 0, 0));
  }

  int N(int lc, int rc, int val) {
    t[cnt] = Node(lc, rc, val);
    return cnt++;
  }

  int ins(int p, int pl, int pr, int x) {
    if (pl > x || pr < x) return p;
    if (pl == pr) return N(0, 0, t[p].val + 1);
    return N(ins(t[p].lc, pl, mid, x), ins(t[p].rc, mid + 1, pr, x), t[p].val + 1);
  }

  int ask(int p1, int p2, int pl, int pr, int k) {
    if (pl == pr) return pl;
    i64 vl = t[t[p2].lc].val - t[t[p1].lc].val;
    if (k <= vl) return ask(t[p1].lc, t[p2].lc, pl, mid, k);
    return ask(t[p1].rc, t[p2].rc, mid + 1, pr, k - vl);
  }

  void add(int x) {
    root.push_back(ins(root.back(), 0, size - 1, x));
  }

  int query(int l, int r, int k) {
    return ask(root[l - 1], root[r], 0, size - 1, k);
  }

#undef mid
};
```

### pb_ds

```cpp
// 平衡树
#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;
template<class T>
using rank_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template<class Key, class T>
using rank_map = tree<Key, T, less<Key>, rb_tree_tag, tree_order_statistics_node_update>;

void example() {
  rank_set<int> t, t2;
  t.insert(8);
  auto it = t.insert(10).first;
  assert(it == t.lower_bound(9));
  assert(t.order_of_key(10) == 1);
  assert(t.order_of_key(11) == 2);
  assert(*t.find_by_order(0) == 8);
  t.join(t2);  // assuming T < T2 or T > T2, merge t2 into t
}

// 优先队列
#include <ext/pb_ds/priority_queue.hpp>
using namespace __gnu_pbds;
template<class T, class Cmp = less<T> >
using pair_heap = __gnu_pbds::priority_queue<T, Cmp>;

void example() {
  pair_heap<int> q1, q2;
  q1.push(5);
  q1.push(10);
  q2.push(1);
  q2.push(7);
  q1.join(q2);
}
```

### Splay

```cpp
// 正常Splay
struct Node {
  int val, size;
  Node *pa, *lc, *rc;
  Node(int val = 0, Node *pa = nullptr) : val(val), size(1), pa(pa), lc(nullptr), rc(nullptr) {}
  Node*& c(bool x) { return x ? lc : rc; }
  bool d() { return pa ? this == pa->lc : 0; }
} pool[MAXN << 2], *tail = pool;

struct Splay {
  Node *root;

  Splay() : root(nullptr) {}

  Node* N(int val, Node *pa) {
    return new (tail++) Node(val, pa);
  }

  void upd(Node *o) {
    o->size = (o->lc ? o->lc->size : 0) + (o->rc ? o->rc->size : 0) + 1;
  }

  void link(Node *x, Node *y, bool d) {
    if (x) x->pa = y;
    if (y) y->c(d) = x;
  }

  void rotate(Node *o) {
    bool dd = o->d();
    Node *x = o->pa, *xx = x->pa, *y = o->c(!dd);
    link(o, xx, x->d());
    link(y, x, dd);
    link(x, o, !dd);
    upd(x);
    upd(o);
  }

  void splay(Node *o) {
    for (Node *x = o->pa; x = o->pa, x; rotate(o)) {
      if (x->pa) rotate(o->d() == x->d() ? x : o);
    }
    root = o;
  }
};
```

### Treap

```cpp
// split_x 左侧元素 < x
// split_k 左侧分割出 k 个元素
namespace tr {
  using uint = unsigned int;

  uint rnd() {
    static uint A = 1 << 16 | 3, B = 33333331, C = 1091;
    return C = A * C + B;
  }

  struct Node {
    uint key;
    int val, size;
    Node *lc, *rc;
    Node(int val = 0) : key(rnd()), val(val), size(1), lc(nullptr), rc(nullptr) {}
  } pool[MAXN << 2], *tail = pool;

  Node* N(int val) {
    return new (tail++) Node(val);
  }

  void upd(Node *o) {
    o->size = (o->lc ? o->lc->size : 0) + (o->rc ? o->rc->size : 0) + 1;
  }

  Node* merge(Node *l, Node *r) {
    if (!l) return r;
    if (!r) return l;
    if (l->key > r->key) {
      l->rc = merge(l->rc, r);
      upd(l);
      return l;
    } else {
      r->lc = merge(l, r->lc);
      upd(r);
      return r;
    }
  }

  void split_x(Node *o, int x, Node*& l, Node*& r) {
    if (!o) { l = r = nullptr; return; }
    if (o->val < x) {
      l = o;
      split_x(o->rc, x, l->rc, r);
      upd(l);
    } else {
      r = o;
      split_x(o->lc, x, l, r->lc);
      upd(r);
    }
  }

  void split_k(Node *o, int k, Node*& l, Node*& r) {
    if (!o) { l = r = nullptr; return; }
    int lsize = o->lc ? o->lc->size : 0;
    if (lsize < k) {
      l = o;
      split_k(o->rc, k - lsize - 1, o->rc, r);
      upd(l);
    } else {
      r = o;
      split_k(o->lc, k, l, o->lc);
      upd(r);
    }
  }
}
```

### 珂朵莉树

```cpp
// for Q assign operation it takes Qlogn time in total
template <class T>
struct interval_set {
  map<pair<int, int>, T> mp;  // {r,l}=val
  interval_set(int l, int r, T v = T()) { mp[{r, l}] = v; }

  // assign a[i]=val for l<=i<=r
  // returns affected ranges before performing this assign operation
  vector<pair<pair<int, int>, T>> assign(int l, int r, T v) {
    auto b = mp.lower_bound({l, 0})->first;
    if (b.second != l) {
      T z = mp[b];
      mp.erase(b);
      mp[{l - 1, b.second}] = z;
      mp[{b.first, l}] = z;
    }

    auto e = mp.lower_bound({r, 0})->first;
    if (e.first != r) {
      T z = mp[e];
      mp.erase(e);
      mp[{e.first, r + 1}] = z;
      mp[{r, e.second}] = z;
    }

    vector<pair<pair<int, int>, T>> ret;
    for (auto it = mp.lower_bound({l, 0}); it != mp.end() && it->first.first <= r; ++it) {
      ret.push_back({{it->first.second, it->first.first}, it->second});
    }
    for (auto it : ret) mp.erase({it.first.second, it.first.first});
    mp[{r, l}] = v;
    return ret;
  }

  T operator[](int i) const { return mp.lower_bound({i, 0})->second; }
};
```

### CDQ 分治

+ 三维偏序（不严格）

```cpp
struct Node {
  int x, y, z, sum, ans;
} p[N], q[N];

void CDQ(int l, int r) {
  if (l == r) return;
  int mid = (l + r) >> 1;
  CDQ(l, mid);
  CDQ(mid + 1, r);
  int i = l, j = mid + 1;
  for (int t = l; t <= r; t++) {
    if (j > r || (i <= mid && p[i].y <= p[j].y)) {
      q[t] = p[i++];
      bit.add(q[t].z, q[t].sum);
    } else {
      q[t] = p[j++];
      q[t].ans += bit.get(q[t].z);
    }
  }
  for (i = l; i <= r; i++) {
    p[i] = q[i];
    bit.set(p[i].z, 0);
  }
}

void go() {
  sort(p, p + n, [](const Node &a, const Node &b) {
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    return a.z < b.z;
  });
  auto eq = [](const Node& a, const Node& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  };
  int k = n;
  for (int i = 0, j = 0; i < n; i++, j++) {
    if (eq(p[i], p[j - 1])) j--, k--;
    else if (i != j) p[j] = p[i];
    p[j].sum++;
  }
  bit.init(m);
  CDQ(0, k - 1);
}
```

### 莫队

```cpp
int unit = max(1, int(n / sqrt(q)));
sort(qry.begin(), qry.end(), [&](auto &a, auto &b) {
  if (a.l / unit != b.l / unit) return a.l < b.l;
  return ((a.l / unit) & 1) ? a.r < b.r : a.r > b.r;
});

// [l, r)
while (l > q.l) add_left(--l);
while (r < q.r) add_right(r++);
while (l < q.l) delete_left(l++);
while (r > q.r) delete_right(--r);
```