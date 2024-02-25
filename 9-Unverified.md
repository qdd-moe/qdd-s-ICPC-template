## Unverified

**Copyright belongs to the original author. Some code has style adjustments. Not guaranteed to be correct.**

### Segment Tree Beats

```cpp
// Nyaan
struct AngelBeats {
  static constexpr i64 INF = numeric_limits<i64>::max() / 2.1;

  struct alignas(32) Node {
    i64 sum = 0, g1 = 0, l1 = 0;
    i64 g2 = -INF, gc = 1, l2 = INF, lc = 1, add = 0;
  };

  vector<Node> v;
  i64 n, log;

  AngelBeats() {}
  AngelBeats(int _n) : AngelBeats(vector<i64>(_n)) {}
  AngelBeats(const vector<i64>& vc) {
    n = 1, log = 0;
    while (n < (int)vc.size()) n <<= 1, log++;
    v.resize(2 * n);
    for (i64 i = 0; i < (int)vc.size(); ++i) {
      v[i + n].sum = v[i + n].g1 = v[i + n].l1 = vc[i];
    }
    for (i64 i = n - 1; i; --i) update(i);
  }

  void range_chmin(int l, int r, i64 x) { inner_apply<1>(l, r, x); }
  void range_chmax(int l, int r, i64 x) { inner_apply<2>(l, r, x); }
  void range_add(int l, int r, i64 x) { inner_apply<3>(l, r, x); }
  void range_update(int l, int r, i64 x) { inner_apply<4>(l, r, x); }
  i64 range_min(int l, int r) { return inner_fold<1>(l, r); }
  i64 range_max(int l, int r) { return inner_fold<2>(l, r); }
  i64 range_sum(int l, int r) { return inner_fold<3>(l, r); }

 private:
  void update(int k) {
    Node& p = v[k];
    Node& l = v[k * 2 + 0];
    Node& r = v[k * 2 + 1];

    p.sum = l.sum + r.sum;

    if (l.g1 == r.g1) {
      p.g1 = l.g1;
      p.g2 = max(l.g2, r.g2);
      p.gc = l.gc + r.gc;
    } else {
      bool f = l.g1 > r.g1;
      p.g1 = f ? l.g1 : r.g1;
      p.gc = f ? l.gc : r.gc;
      p.g2 = max(f ? r.g1 : l.g1, f ? l.g2 : r.g2);
    }

    if (l.l1 == r.l1) {
      p.l1 = l.l1;
      p.l2 = min(l.l2, r.l2);
      p.lc = l.lc + r.lc;
    } else {
      bool f = l.l1 < r.l1;
      p.l1 = f ? l.l1 : r.l1;
      p.lc = f ? l.lc : r.lc;
      p.l2 = min(f ? r.l1 : l.l1, f ? l.l2 : r.l2);
    }
  }

  void push_add(int k, i64 x) {
    Node& p = v[k];
    p.sum += x << (log + __builtin_clz(k) - 31);
    p.g1 += x;
    p.l1 += x;
    if (p.g2 != -INF) p.g2 += x;
    if (p.l2 != INF) p.l2 += x;
    p.add += x;
  }
  void push_min(int k, i64 x) {
    Node& p = v[k];
    p.sum += (x - p.g1) * p.gc;
    if (p.l1 == p.g1) p.l1 = x;
    if (p.l2 == p.g1) p.l2 = x;
    p.g1 = x;
  }
  void push_max(int k, i64 x) {
    Node& p = v[k];
    p.sum += (x - p.l1) * p.lc;
    if (p.g1 == p.l1) p.g1 = x;
    if (p.g2 == p.l1) p.g2 = x;
    p.l1 = x;
  }
  void push(int k) {
    Node& p = v[k];
    if (p.add != 0) {
      push_add(k * 2 + 0, p.add);
      push_add(k * 2 + 1, p.add);
      p.add = 0;
    }
    if (p.g1 < v[k * 2 + 0].g1) push_min(k * 2 + 0, p.g1);
    if (p.l1 > v[k * 2 + 0].l1) push_max(k * 2 + 0, p.l1);

    if (p.g1 < v[k * 2 + 1].g1) push_min(k * 2 + 1, p.g1);
    if (p.l1 > v[k * 2 + 1].l1) push_max(k * 2 + 1, p.l1);
  }

  void subtree_chmin(int k, i64 x) {
    if (v[k].g1 <= x) return;
    if (v[k].g2 < x) {
      push_min(k, x);
      return;
    }
    push(k);
    subtree_chmin(k * 2 + 0, x);
    subtree_chmin(k * 2 + 1, x);
    update(k);
  }

  void subtree_chmax(int k, i64 x) {
    if (x <= v[k].l1) return;
    if (x < v[k].l2) {
      push_max(k, x);
      return;
    }
    push(k);
    subtree_chmax(k * 2 + 0, x);
    subtree_chmax(k * 2 + 1, x);
    update(k);
  }

  template <int cmd>
  inline void _apply(int k, i64 x) {
    if constexpr (cmd == 1) subtree_chmin(k, x);
    if constexpr (cmd == 2) subtree_chmax(k, x);
    if constexpr (cmd == 3) push_add(k, x);
    if constexpr (cmd == 4) subtree_chmin(k, x), subtree_chmax(k, x);
  }

  template <int cmd>
  void inner_apply(int l, int r, i64 x) {
    if (l == r) return;
    l += n, r += n;
    for (int i = log; i >= 1; i--) {
      if (((l >> i) << i) != l) push(l >> i);
      if (((r >> i) << i) != r) push((r - 1) >> i);
    }
    {
      int l2 = l, r2 = r;
      while (l < r) {
        if (l & 1) _apply<cmd>(l++, x);
        if (r & 1) _apply<cmd>(--r, x);
        l >>= 1;
        r >>= 1;
      }
      l = l2;
      r = r2;
    }
    for (int i = 1; i <= log; i++) {
      if (((l >> i) << i) != l) update(l >> i);
      if (((r >> i) << i) != r) update((r - 1) >> i);
    }
  }

  template <int cmd>
  inline i64 e() {
    if constexpr (cmd == 1) return INF;
    if constexpr (cmd == 2) return -INF;
    return 0;
  }

  template <int cmd>
  inline void op(i64& a, const Node& b) {
    if constexpr (cmd == 1) a = min(a, b.l1);
    if constexpr (cmd == 2) a = max(a, b.g1);
    if constexpr (cmd == 3) a += b.sum;
  }

  template <int cmd>
  i64 inner_fold(int l, int r) {
    if (l == r) return e<cmd>();
    l += n, r += n;
    for (int i = log; i >= 1; i--) {
      if (((l >> i) << i) != l) push(l >> i);
      if (((r >> i) << i) != r) push((r - 1) >> i);
    }
    i64 lx = e<cmd>(), rx = e<cmd>();
    while (l < r) {
      if (l & 1) op<cmd>(lx, v[l++]);
      if (r & 1) op<cmd>(rx, v[--r]);
      l >>= 1;
      r >>= 1;
    }
    if constexpr (cmd == 1) lx = min(lx, rx);
    if constexpr (cmd == 2) lx = max(lx, rx);
    if constexpr (cmd == 3) lx += rx;
    return lx;
  }
};
```

### Josephus Problem

```cpp
// n people, count from 1 to m, asking for the number of the last person remaining
// Formula: f(n,m)=(f(n−1,m)+m)%n, f(0,m)=0;
// O(n)
i64 calc(int n, i64 m) {
    i64 p = 0;
    for (int i = 2; i <= n; i++) {
        p = (p + m) % i;
    }
    return p + 1;
}

// n people, count from 1 to m, asking for the number of the k-th person eliminated
// Formula: f(n,k)=(f(n−1,k−1)+m−1)%n+1
// f(n−k+1,1)=m%(n−k+1)
// if (f==0) f=n−k+1
// O(k)
i64 cal1(i64 n, i64 m, i64 k) {  // (k == n) equal(calc)
    i64 p = m % (n - k + 1);
    if (p == 0) p = n - k + 1;
    for (i64 i = 2; i <= k; i++) {
        p = (p + m - 1) % (n - k + i) + 1;
    }
    return p;
}

// n people, count from 1 to m, asking for the number of the k-th person eliminated
// O(m*log(m))
i64 cal2(i64 n, i64 m, i64 k) {
    if (m == 1)
        return k;
    else {
        i64 a = n - k + 1, b = 1;
        i64 c = m % a, x = 0;
        if (c == 0) c = a;
        while (b + x <= k) {
            a += x, b += x, c += m * x;
            c %= a;
            if (c == 0) c = a;
            x = (a - c) / (m - 1) + 1;
        }
        c += (k - b) * m;
        c %= n;
        if (c == 0) c = n;
        return c;
    }
}

// n people, count from 1 to m, asking for the number of the person with index k eliminated
// O(n)
i64 n, k;  // n <= 4e7, number of queries <= 100, index range [0,n-1]
i64 dieInXturn(int n, int k, int x) {  // n people, count k, the X-th person dies with index X
    i64 tmp = 0;
    while (n) {
        x = (x + n) % n;
        if (k > n) x += (k - x - 1 + n - 1) / n * n;
        if ((x + 1) % k == 0) {
            tmp += (x + 1) / k;
            break;
        } else {
            if (k > n) {
                tmp += x / k;
                i64 ttmp = x;
                x = x - (x / n + 1) * (x / k) + (x + n) / n * n - k;
                n -= ttmp / k;
            } else {
                tmp += n / k;
                x = x - x / k;
                x += n - n / k * k;
                n -= n / k;
            }
        }
    }
    return tmp;
}
```

### Lexicographically Smallest Solution of 2-SAT

```cpp
const int N = 1e5 + 10;
struct TwoSatBF {  // Brute-force to find the lexicographically smallest solution
  int n;
  vector<int> G[N << 1];
  bool slt[N << 1];
  // Even points: false Odd points: true So x^1 is the opposite side
  void init(int _n) {
    n = _n;
    for (int i = 0; i < (n << 1); ++i) {
      G[i].clear();
      slt[i] = false;
    }
  }
  void addLimit(int x, int y) {
    // If x is chosen, y must be chosen as well, depending on the situation
    G[x].push_back(y);
    G[y ^ 1].push_back(x ^ 1);
  }
  stack<int> st;
  void clearst() {
    while (st.size()) st.pop();
  }
  bool dfs(int u) {
    if (slt[u ^ 1]) {
      return false;
    } else if (slt[u]) {
      return true;
    }
    slt[u] = true;
    st.push(u);
    for (auto v : G[u]) {
      if (!dfs(v)) {
        return false;
      }
    }
    return true;
  }
  bool solve() {
    for (int u = 0; u < (n << 1); u += 2) {
      if (!slt[u] && !slt[u ^ 1]) {
        clearst();
        if (!dfs(u)) {
          clearst();
          if (!dfs(u ^ 1)) {
            return false;
          }
        }
      }
    }
    return true;
  }
};
```

### Weighted Bipartite Matching (KM Algorithm)

```cpp
// ECNU
namespace R {
    int n;
    int w[N][N], kx[N], ky[N], py[N], vy[N], slk[N], pre[N];

    i64 go() {
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++)
                kx[i] = max(kx[i], w[i][j]);
        for (int i = 1; i <= n; i++) {
            fill(vy, vy + n + 1, 0);
            fill(slk, slk + n + 1, INF);
            fill(pre, pre + n + 1, 0);
            int k = 0, p = -1;
            for (py[k = 0] = i; py[k]; k = p) {
                int d = INF;
                vy[k] = 1;
                int x = py[k];
                for (int j = 1; j <= n; j++) {
                    if (!vy[j]) {
                        int t = kx[x] + ky[j] - w[x][j];
                        if (t < slk[j]) { slk[j] = t; pre[j] = k; }
                        if (slk[j] < d) { d = slk[j]; p = j; }
                    }
                }
                for (int j = 0; j <= n; j++) {
                    if (vy[j]) { kx[py[j]] -= d; ky[j] += d; }
                    else slk[j] -= d;
                }
            }
            for (; k; k = pre[k]) py[k] = py[pre[k]];
        }
        i64 ans = 0;
        for (int i = 1; i <= n; i++) ans += kx[i] + ky[i];
        return ans;
    }
}
```

### HLPP

```cpp
struct HLPP {
    struct Edge {
        int v, rev;
        i64 cap;
    };
    int n, sp, tp, lim, ht, lcnt;
    i64 exf[N];
    vector<Edge> G[N];
    vector<int> hq[N], gap[N], h, sum;
    void init(int nn, int s, int t) {
        sp = s, tp = t, n = nn, lim = n + 1, ht = lcnt = 0;
        for (int i = 1; i <= n; ++i) G[i].clear(), exf[i] = 0;
    }
    void add_edge(int u, int v, i64 cap) {
        G[u].push_back({v, int(G[v].size()), cap});
        G[v].push_back({u, int(G[u].size()) - 1, 0});
    }
    void update(int u, int nh) {
        ++lcnt;
        if (h[u] != lim) --sum[h[u]];
        h[u] = nh;
        if (nh == lim) return;
        ++sum[ht = nh];
        gap[nh].push_back(u);
        if (exf[u] > 0) hq[nh].push_back(u);
    }
    void relabel() {
        queue<int> q;
        for (int i = 0; i <= lim; ++i) hq[i].clear(), gap[i].clear();
        h.assign(lim, lim), sum.assign(lim, 0), q.push(tp);
        lcnt = ht = h[tp] = 0;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (Edge& e : G[u])
                if (h[e.v] == lim && G[e.v][e.rev].cap) update(e.v, h[u] + 1), q.push(e.v);
            ht = h[u];
        }
    }
    void push(int u, Edge& e) {
        if (!exf[e.v]) hq[h[e.v]].push_back(e.v);
        i64 df = min(exf[u], e.cap);
        e.cap -= df, G[e.v][e.rev].cap += df;
        exf[u] -= df, exf[e.v] += df;
    }
    void discharge(int u) {
        int nh = lim;
        if (h[u] == lim) return;
        for (Edge& e : G[u]) {
            if (!e.cap) continue;
            if (h[u] == h[e.v] + 1) {
                push(u, e);
                if (exf[u] <= 0) return;
            } else if (nh > h[e.v] + 1)
                nh = h[e.v] + 1;
        }
        if (sum[h[u]] > 1)
            update(u, nh);
        else {
            for (; ht >= h[u]; gap[ht--].clear())
                for (int& i : gap[ht]) update(i, lim);
        }
    }
    i64 hlpp() {
        exf[sp] = INF, exf[tp] = -INF, relabel();
        for (Edge& e : G[sp]) push(sp, e);
        for (; ~ht; --ht) {
            while (!hq[ht].empty()) {
                int u = hq[ht].back();
                hq[ht].pop_back();
                discharge(u);
                if (lcnt > (n << 2)) relabel();
            }
        }
        return exf[tp] + INF;
    }
};
```

### Flow with Lower Bounds

```cpp
const int INF = 0x3f3f3f3f;

struct edge {
    int to, cap, rev;
};

const int N = 60003;
const int M = 400003;

struct graph {
    int n, m;
    edge w[M];
    int fr[M];
    int num[N], cur[N], first[N];
    edge e[M];

    void init(int n) {
        this->n = n;
        m = 0;
    }

    void add_edge(int from, int to, int cap) {
        w[++m] = (edge){to, cap};
        num[from]++, fr[m] = from;
        w[++m] = (edge){from, 0};
        num[to]++, fr[m] = to;
    }

    void prepare() {
        first[1] = 1;
        for (int i = 2; i <= n; i++) first[i] = first[i - 1] + num[i - 1];
        for (int i = 1; i < n; i++) num[i] = first[i + 1] - 1;
        num[n] = m;
        for (int i = 1; i <= m; i++) {
            e[first[fr[i]] + (cur[fr[i]]++)] = w[i];

            if (!(i % 2)) {
                e[first[fr[i]] + cur[fr[i]] - 1].rev =
                    first[w[i].to] + cur[w[i].to] - 1;
                e[first[w[i].to] + cur[w[i].to] - 1].rev =
                    first[fr[i]] + cur[fr[i]] - 1;
            }
        }
    }

    int q[N];
    int dist[N];
    int t;

    bool bfs(int s) {
        int l = 1, r = 1;
        q[1] = s;
        memset(dist, -1, (n + 1) * 4);
        dist[s] = 0;
        while (l <= r) {
            int u = q[l++];
            for (int i = first[u]; i <= num[u]; i++) {
                int v = e[i].to;
                if ((dist[v] != -1) || (!e[i].cap)) continue;
                dist[v] = dist[u] + 1;
                if (v == t) return true;
                q[++r] = v;
            }
        }
        return dist[t] != -1;
    }

    int dfs(int u, int flow) {
        if (u == t) return flow;
        for (int& i = cur[u]; i <= num[u]; i++) {
            int v = e[i].to;
            if (!e[i].cap || dist[v] != dist[u] + 1) continue;
            int t = dfs(v, min(flow, e[i].cap));
            if (t) {
                e[i].cap -= t;
                e[e[i].rev].cap += t;
                return t;
            }
        }
        return 0;
    }

    i64 dinic(int s, int t) {
        i64 ans = 0;
        this->t = t;
        while (bfs(s)) {
            int flow;
            for (int i = 1; i <= n; i++) cur[i] = first[i];
            while (flow = dfs(s, INF)) ans += (i64)flow;
        }
        return ans;
    }
};

struct graph_bounds {
    int in[N];
    int S, T, sum, cur;
    graph g;
    int n;

    void init(int n) {
        this->n = n;
        S = n + 1;
        T = n + 2;
        sum = 0;
        g.init(n + 2);
    }

    void add_edge(int from, int to, int low, int up) {
        g.add_edge(from, to, up - low);
        in[to] += low;
        in[from] -= low;
    }

    void build() {
        for (int i = 1; i <= n; i++)
            if (in[i] > 0)
                g.add_edge(S, i, in[i]), sum += in[i];
            else if (in[i])
                g.add_edge(i, T, -in[i]);
        g.prepare();
    }

    bool canflow() {
        build();
        int flow = g.dinic(S, T);
        return flow >= sum;
    }

    bool canflow(int s, int t) {
        g.add_edge(t, s, INF);
        build();
        for (int i = 1; i <= g.m; i++) {
            edge& e = g.e[i];
            if (e.to == s && e.cap == INF) {
                cur = i;
                break;
            }
        }
        int flow = g.dinic(S, T);
        return flow >= sum;
    }

    int maxflow(int s, int t) {
        if (!canflow(s, t)) return -1;
        return g.dinic(s, t);
    }

    int minflow(int s, int t) {
        if (!canflow(s, t)) return -1;
        edge& e = g.e[cur];
        int flow = INF - e.cap;
        e.cap = g.e[e.rev].cap = 0;
        return flow - g.dinic(t, s);
    }
} g;

void solve() {
    int n = read(), m = read(), s = read(), t = read();
    g.init(n);
    while (m--) {
        int u = read(), v = read(), low = read(), up = read();
        g.add_edge(u, v, low, up);
    }
}
```

### Heavy-Light Decomposition

```cpp
// jiangly
struct HLD {
    int n;
    vector<int> siz, top, dep, parent, in, out, seq;
    vector<vector<int>> adj;
    int cur;

    explicit HLD(int n) : n(n), siz(n), top(n), dep(n), parent(n), in(n), out(n), seq(n), adj(n), cur(0) {}

    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    void work(int root = 0) {
        top[root] = root;
        dep[root] = 0;
        parent[root] = -1;
        dfs1(root);
        dfs2(root);
    }

    void dfs1(int u) {
        if (parent[u] != -1) {
            adj[u].erase(find(adj[u].begin(), adj[u].end(), parent[u]));
        }
        siz[u] = 1;
        for (auto& v : adj[u]) {
            parent[v] = u;
            dep[v] = dep[u] + 1;
            dfs1(v);
            siz[u] += siz[v];
            if (siz[v] > siz[adj[u][0]]) {
                swap(v, adj[u][0]);
            }
        }
    }

    void dfs2(int u) {
        in[u] = cur++;
        seq[in[u]] = u;
        for (auto v : adj[u]) {
            top[v] = v == adj[u][0] ? top[u] : v;
            dfs2(v);
        }
        out[u] = cur;
    }

    int lca(int u, int v) {
        while (top[u] != top[v]) {
            if (dep[top[u]] > dep[top[v]]) {
                u = parent[top[u]];
            } else {
                v = parent[top[v]];
            }
        }
        return dep[u] < dep[v] ? u : v;
    }

    int dist(int u, int v) { return dep[u] + dep[v] - 2 * dep[lca(u, v)]; }

    int jump(int u, int k) {
        if (dep[u] < k) {
            return -1;
        }
        int d = dep[u] - k;
        while (dep[top[u]] > d) {
            u = parent[top[u]];
        }
        return seq[in[u] - dep[u] + d];
    }

    bool isAncester(int u, int v) { return in[u] <= in[v] && in[v] < out[u]; }

    int rootedParent(int u, int v) {
        swap(u, v);
        if (u == v) {
            return u;
        }
        if (!isAncester(u, v)) {
            return parent[u];
        }
        auto it = upper_bound(adj[u].begin(), adj[u].end(), v, [&](int x, int y) { return in[x] < in[y]; }) - 1;
        return *it;
    }

    int rootedSize(int u, int v) {
        if (u == v) {
            return n;
        }
        if (!isAncester(v, u)) {
            return siz[v];
        }
        return n - siz[rootedParent(u, v)];
    }

    int rootedLca(int a, int b, int c) { return lca(a, b) ^ lca(b, c) ^ lca(c, a); }
};
```

### Link-Cut Tree

```cpp
// Chestnut
const int N = 50005;

#define lc son[x][0]
#define rc son[x][1]

struct Splay {
    int fa[N], son[N][2];
    int st[N];
    bool rev[N];
    inline int which(int x) {
        for (int i = 0; i < 2; i++)
            if (son[fa[x]][i] == x) return i;
        return -1;
    }

    inline void pushdown(int x) {
        if (rev[x]) {
            rev[x] ^= 1;
            rev[lc] ^= 1;
            rev[rc] ^= 1;
            swap(lc, rc);
        }
    }
    
    inline void rotate(int x) {
        int f = fa[x], w = which(x) ^ 1, c = son[x][w];
        fa[x] = fa[f];
        if (which(f) != -1) son[fa[f]][which(f)] = x;
        fa[c] = f;
        son[f][w ^ 1] = c;
        fa[f] = x;
        son[x][w] = f;
    }

    inline void splay(int x) {
        int top = 0;
        st[++top] = x;
        for (int i = x; which(i) != -1; i = fa[i]) {
            st[++top] = fa[i];
        }
        for (int i = top; i; i--) pushdown(st[i]);
        while (which(x) != -1) {
            int f = fa[x];
            if (which(f) != -1) {
                if (which(x) ^ which(f)) rotate(x);
                else rotate(f);
            }
            rotate(x);
        }
    }

    void access(int x) {
        int t = 0;
        while (x) {
            splay(x);
            rc = t;
            t = x;
            x = fa[x];
        }
    }

    void rever(int x) {
        access(x);
        splay(x);
        rev[x] ^= 1;
    }

    void link(int x, int y) {
        rever(x);
        fa[x] = y;
        splay(x);
    }

    void cut(int x, int y) {
        rever(x);
        access(y);
        splay(y);
        son[y][0] = fa[x] = 0;
    }

    int find(int x) {
        access(x);
        splay(x);
        int y = x;
        while (son[y][0]) y = son[y][0];
        return y;
    }
} T;

int n, m;

int main() {
    char ch[10];
    int x, y;
    scanf("%d%d", &n, &m);
    for (int i = 1; i <= m; i++) {
        scanf("%s", ch);
        scanf("%d%d", &x, &y);
        if (ch[0] == 'C') T.link(x, y);
        else if (ch[0] == 'D') T.cut(x, y);
        else {
            if (T.find(x) == T.find(y)) printf("Yes\n");
            else printf("No\n");
        }
    }
}
```

### NTT for NTT-unfriendly Modulus

```cpp
// memset0
const int N = 4e5 + 10, G = 3, P[3] = {469762049, 998244353, 1004535809};
int n1, n2, k, n, p, p1, p2, M2;
int a[N], b[N], f[3][N], g[N], rev[N], ans[N];

void ntt(int *a, int g, int p) {
    for (int i = 0; i < n; i++) if (i < rev[i]) swap(a[i], a[rev[i]]);
    for (int len = 1; len < n; len <<= 1) {
        int wn = qk(g, (p - 1) / (len << 1), p);
        for (int i = 0; i < n; i += (len << 1)) {
            int w = 1;
            for (int j = 0; j < len; j++, w = (i64)w * wn % p) {
                int x = a[i + j], y = (i64)w * a[i + j + len] % p;
                a[i + j] = (x + y) % p, a[i + j + len] = (x - y + p) % p;
            }
        }
    }
}

int merge(int a1, int a2, int A2) {
    i64 M1 = (i64)p1 * p2;
    i64 A1 = ((i64)inv(p2, p1) * a1 % p1 * p2 + (i64)inv(p1, p2) * a2 % p2 * p1) % M1;
    i64 K = ((A2 - A1) % M2 + M2) % M2 * inv(M1 % M2, M2) % M2;
    int ans = (A1 + M1 % p * K) % p;
    return ans;
}

void go() {
    read(n1), read(n2), read(p);
    p1 = P[0], p2 = P[1], M2 = P[2];
    for (int i = 0; i <= n1; i++) read(a[i]);
    for (int i = 0; i <= n2; i++) read(b[i]);
    n = 1; while (n <= (n1 + n2)) n <<= 1, ++k;
    for (int i = 0; i < n; i++) {
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (k - 1));
    }
    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < n; i++) f[k][i] = a[i] % P[k];
        for (int i = 0; i < n; i++) g[i] = b[i] % P[k];
        ntt(f[k], G, P[k]), ntt(g, G, P[k]);
        for (int i = 0; i < n; i++) f[k][i] = (i64)f[k][i] * g[i] % P[k];
        ntt(f[k], inv(G, P[k]), P[k]);
        for (int i = 0; i < n; i++) f[k][i] = (i64)f[k][i] * inv(n, P[k]) % P[k];
    }
    for (int i = 0; i <= n1 + n2; i++) ans[i] = merge(f[0][i], f[1][i], f[2][i]);
}
```

### Geometry

```cpp
// Great-circle distance on the sphere
// Voleking
ld Dist(ld la1, ld lo1, ld la2, ld lo2, ld R) {
    la1 *= PI / 180, lo1 *= PI / 180, la2 *= PI / 180, lo2 *= PI / 180;
    ld x1 = cos(la1) * sin(lo1), y1 = cos(la1) * cos(lo1), z1 = sin(la1); 
    ld x2 = cos(la2) * sin(lo2), y2 = cos(la2) * cos(lo2), z1 = sin(la2); 
    return R * acos(x1 * x2 + y1 * y2 + z1 * z2);
}

// jiry_2
int cmp(ld k1, ld k2) {
    return sgn(k1 - k2);
}
V proj(V k1, V k2, V q) { // Projection of point q onto line k1,k2
    V k = k2 - k1;
    return k1 + k * (dot(q - k1, k) / k.abs2());
}
V reflect(V k1, V k2, V q) {
    return proj(k1, k2, q) * 2 - q;
}
int clockwise(V k1, V k2, V k3) { // k1 k2 k3 counterclockwise 1 clockwise -1 otherwise 0  
    return sgn(det(k2 - k1, k3 - k1));
}
int checkLL(V k1, V k2, V k3, V k4) { // Check the intersection point of line (L) and segment (S) k1,k2 and k3,k4
    return cmp(det(k3 - k1, k4 - k1), det(k3 - k2, k4 - k2)) != 0;
}
V getLL(V k1, V k2, V k3, V k4) {
    ld w1 = det(k1 - k3, k4 - k3), w2 = det(k4 - k3, k2 - k3);
    return (k1 * w2 + k2 * w1) / (w1 + w2);
}
vector<line> getHL(vector<line>& L) { // Get the half-plane intersection, the half-plane is counterclockwise, and the output is counterclockwise
    sort(L.begin(), L.end());
    deque<line> q;
    for (int i = 0; i < (int) L.size(); i++) {
        if (i && sameDir(L[i], L[i - 1])) continue;
        while (q.size() > 1 && !checkpos(q[q.size() - 2], q[q.size() - 1], L[i])) q.pop_back();
        while (q.size() > 1 && !checkpos(q[1], q[0], L[i])) q.pop_front();
        q.push_back(L[i]);
    }
    while (q.size() > 2 && !checkpos(q[q.size() - 2], q[q.size() - 1], q[0])) q.pop_back();
    while (q.size() > 2 && !checkpos(q[1], q[0], q[q.size() - 1])) q.pop_front();
    vector<line> ans;
    for (int i = 0; i < q.size(); i++) ans.push_back(q[i]);
    return ans;
}
int checkposCC(circle k1, circle k2) { // Return the number of common tangent lines of two circles
    if (cmp(k1.r, k2.r) == -1) swap(k1, k2);
    ld dis = k1.o.dis(k2.o);
    int w1 = cmp(dis, k1.r + k2.r), w2 = cmp(dis, k1.r - k2.r);
    if (w1 > 0) return 4;
    else if (w1 == 0) return 3;
    else if (w2 > 0) return 2;
    else if (w2 == 0) return 1;
    else return 0;
}
vector<V> getCL(circle k1, V k2, V k3) { // Given k2->k3 direction, give out p, the two tangent points
    V k = proj(k2, k3, k1.o);
    ld d = k1.r * k1.r - (k - k1.o).abs2();
    if (sgn(d) == -1) return {};
    V del = (k3 - k2).unit() * sqrt(max((ld) 0.0, d));
    return {k - del, k + del};
}
vector<line> TangentoutCC(circle k1, circle k2) {
    int pd = checkposCC(k1, k2);
    if (pd == 0) return {};
    if (pd == 1) {
        V k = getCC(k1, k2)[0];
        return { (line){k, k} };
    }
    if (cmp(k1.r, k2.r) == 0) {
        V del = (k2.o - k1.o).unit().turn90().getdel();
        return {
            (line){k1.o - del * k1.r, k2.o - del * k2.r},
            (line){k1.o + del * k1.r, k2.o + del * k2.r}
        };
    } else {
        V p = (k2.o * k1.r - k1.o * k2.r) / (k1.r - k2.r);
        vector<V> A = TangentCP(k1, p), B = TangentCP(k2, p);
        vector<line> ans;
        for (int i = 0; i < A.size(); i++) ans.push_back((line){A[i], B[i]});
        return ans;
    }
}
vector<line> TangentinCC(circle k1, circle k2) {
    int pd = checkposCC(k1, k2);
    if (pd <= 2) return {};
    if (pd == 3) {
        V k = getCC(k1, k2)[0];
        return { (line){k, k} };
    }
    V p = (k2.o * k1.r + k1.o * k2.r) / (k1.r + k2.r);
    vector<V> A = TangentCP(k1, p), B = TangentCP(k2, p);
    vector<line> ans;
    for (int i = 0; i < A.size(); i++) ans.push_back((line){A[i], B[i]});
    return ans;
}
vector<line> TangentCC(circle k1, circle k2) {
    int flag = 0;
    if (k1.r < k2.r) swap(k1, k2), flag = 1;
    vector<line> A = TangentoutCC(k1, k2), B = TangentinCC(k1, k2);
    for (line k: B) A.push_back(k);
    if (flag) for (line& k: A) swap(k[0], k[1]);
    return A;
}
ld convexDiameter(vector<V> A) {
    int now = 0, n = A.size();
    ld ans = 0;
    for (int i = 0; i < A.size(); i++) {
        now = max(now, i);
        while (1) {
            ld k1 = A[i].dis(A[now % n]), k2 = A[i].dis(A[(now + 1) % n]);
            ans = max(ans, max(k1, k2));
            if (k2 > k1) now++;
            else break;
        }
    }
    return ans;
}
vector<V> convexcut(vector<V> A, V k1, V k2) { // Keep points k1,k2,p counterclockwise
    int n = A.size();
    A.push_back(A[0]);
    vector<V> ans;
    for (int i = 0; i < n; i++) {
        int w1 = clockwise(k1, k2, A[i]), w2 = clockwise(k1, k2, A[i + 1]);
        if (w1 >= 0) ans.push_back(A[i]);
        if (w1 * w2 < 0) ans.push_back(getLL(k1, k2, A[i], A[i + 1]));
    }
    return ans;
}
```

### References

[Nyaan's Library](https://nyaannyaan.github.io/library/)
[F0RE1GNERS](https://github.com/F0RE1GNERS/template)
[YouKn0wWho](https://github.com/ShahjalalShohag/code-library)
[OI Wiki](https://oi-wiki.org/)
[cp-algorithms](https://cp-algorithms.com/)
