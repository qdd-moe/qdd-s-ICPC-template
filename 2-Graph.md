## Graph Theory

### Adjacency List

```cpp
int ecnt, mp[N];

struct Edge {
  int to, nxt;
  Edge(int to = 0, int nxt = 0) : to(to), nxt(nxt) {}
} es[M];

void mp_init() {
  memset(mp, -1, (n + 2) * sizeof(int));
  ecnt = 0;
}

void mp_link(int u, int v) {
  es[ecnt] = Edge(v, mp[u]);
  mp[u] = ecnt++;
}

for (int i = mp[u]; i != -1; i = es[i].nxt)
```

### Shortest Path

+ Dijkstra

```cpp
struct Dijkstra {
  struct Edge {
    int to, val;
    Edge(int to = 0, int val = 0) : to(to), val(val) {}
  };

  int n;
  vector<vector<Edge>> g;

  Dijkstra(int n) : n(n), g(n) {}

  void add_edge(int u, int v, int val) { g[u].emplace_back(v, val); }

  vector<i64> solve(int s) {
    using pii = pair<i64, int>;
    vector<i64> dis(n, 1LL << 60);
    priority_queue<pii, vector<pii>, greater<pii>> q;
    dis[s] = 0;
    q.emplace(0, s);
    while (!q.empty()) {
      pii p = q.top();
      q.pop();
      int u = p.second;
      if (dis[u] < p.first) continue;
      for (Edge& e : g[u]) {
        int v = e.to;
        if (dis[v] > dis[u] + e.val) {
          dis[v] = dis[u] + e.val;
          q.emplace(dis[v], v);
        }
      }
    }
    return dis;
  }
};
```

+ Bellman-Ford

```cpp
struct SPFA {
  struct Edge {
    int to, val;
    Edge(int to = 0, int val = 0) : to(to), val(val) {}
  };

  int n;
  vector<vector<Edge>> g;

  SPFA(int n) : n(n), g(n) {}

  void add_edge(int u, int v, int val) { g[u].emplace_back(v, val); }

  vector<i64> solve(int s) {
    queue<int> q;
    vector<i64> dis(n, 1LL << 60);
    vector<bool> in(n, 0);
    q.push(s);
    dis[s] = 0;
    in[s] = 1;
    while (!q.empty()) {
      int u = q.front();
      q.pop();
      in[u] = 0;
      for (Edge& e : g[u]) {
        int v = e.to;
        if (dis[v] > dis[u] + e.val) {
          dis[v] = dis[u] + e.val;
          if (!in[v]) {
            in[v] = 1;
            q.push(v);
          }
        }
      }
    }
    return dis;
  }
};
```

+ Floyd, with Shortest Cycle

```cpp
// Note: INF should not exceed 1/3 LLONG_MAX
for (int k = 0; k < n; k++) {
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < i; j++) {
      ans = min(ans, g[i][k] + g[k][j] + dis[i][j]);
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      dis[i][j] = min(dis[i][j], dis[i][k] + dis[k][j]);
    }
  }
}
```

### Topological Sorting

```cpp
// Use max/min heap for lexicographically largest/smallest
int n, deg[N], dis[N];
vector<int> g[N];

bool topo(vector<int>& ans) {
  queue<int> q;
  for (int i = 1; i <= n; i++) {
    if (deg[i] == 0) {
      q.push(i);
      dis[i] = 1;
    }
  }
  while (!q.empty()) {
    int u = q.front();
    q.pop();
    ans.push_back(u);
    for (int v : g[u]) {
      deg[v]--;
      dis[v] = max(dis[v], dis[u] + 1);
      if (deg[v] == 0) q.push(v);
    }
  }
  return ans.size() == n;
}
```

### Minimum Spanning Tree

```cpp
// Prerequisite: Disjoint Set Union
struct Edge {
  int u, v, w;
  Edge(int u = 0, int v = 0, int w = 0) : u(u), v(v), w(w) {}
};

i64 kruskal(vector<Edge>& es, int n) {
  sort(es.begin(), es.end(), [](Edge& x, Edge& y) { return x.w < y.w; });
  dsu d(n + 1);
  i64 ans = 0;
  for (Edge& e : es) {
    if (d.merge(e.u, e.v)) {
      ans += e.w;
    }
  }
  return ans;
}
```

### Lowest Common Ancestor

```cpp
// dfs(1, 0) or dfs(0, n), don't use dfs(0, -1)
const int LOG = 22;  // 22 = ((int)log2(N) + 1)
int dep[N], up[N][LOG];

void dfs(int u, int pa) {
  dep[u] = dep[pa] + 1;
  up[u][0] = pa;
  for (int i = 1; i < LOG; i++) {
    up[u][i] = up[up[u][i - 1]][i - 1];
  }
  for (int v : g[u]) {
    if (v != pa) dfs(v, u);
  }
}

int lca(int u, int v) {
  if (dep[u] > dep[v]) swap(u, v);
  int t = dep[v] - dep[u];
  for (int i = 0; i < LOG; i++) {
    if ((t >> i) & 1) v = up[v][i];
  }
  if (u == v) return u;
  for (int i = LOG - 1; i >= 0; i--) {
    if (up[u][i] != up[v][i]) {
      u = up[u][i];
      v = up[v][i];
    }
  }
  return up[u][0];
}
```

### Network Flow

+ Max Flow

```cpp
const int INF = 0x7fffffff;

struct Dinic {
  struct Edge {
    int to, cap;
    Edge(int to, int cap) : to(to), cap(cap) {}
  };

  int n, s, t;
  vector<Edge> es;
  vector<vector<int>> g;
  vector<int> dis, cur;

  Dinic(int n, int s, int t) : n(n), s(s), t(t), g(n + 1), dis(n + 1), cur(n + 1) {}

  void add_edge(int u, int v, int cap) {
    g[u].push_back(es.size());
    es.emplace_back(v, cap);
    g[v].push_back(es.size());
    es.emplace_back(u, 0);
  }

  bool bfs() {
    dis.assign(n + 1, 0);
    queue<int> q;
    q.push(s);
    dis[s] = 1;
    while (!q.empty()) {
      int u = q.front();
      q.pop();
      for (int i : g[u]) {
        Edge& e = es[i];
        if (!dis[e.to] && e.cap > 0) {
          dis[e.to] = dis[u] + 1;
          q.push(e.to);
        }
      }
    }
    return dis[t];
  }

  int dfs(int u, int cap) {
    if (u == t || cap == 0) return cap;
    int tmp = cap;
    for (int& i = cur[u]; i < (int)g[u].size(); i++) {
      Edge& e = es[g[u][i]];
      if (dis[e.to] == dis[u] + 1) {
        int f = dfs(e.to, min(cap, e.cap));
        e.cap -= f;
        es[g[u][i] ^ 1].cap += f;
        cap -= f;
        if (cap == 0) break;
      }
    }
    return tmp - cap;
  }

  i64 solve() {
    i64 flow = 0;
    while (bfs()) {
      cur.assign(n + 1, 0);
      flow += dfs(s, INF);
    }
    return flow;
  }
};
```

+ Minimum Cost Flow

```cpp
const i64 INF = 1e15;

struct MCMF {
  struct Edge {
    int from, to;
    i64 cap, cost;
    Edge(int from, int to, i64 cap, i64 cost) : from(from), to(to), cap(cap), cost(cost) {}
  };

  int n, s, t;
  i64 flow, cost;
  vector<Edge> es;
  vector<vector<int>> g;
  vector<i64> d, a;  // dis, add, prev
  vector<int> p, in;

  MCMF(int n, int s, int t) : n(n), s(s), t(t), flow(0), cost(0), g(n + 1), p(n + 1), a(n + 1) {}

  void add_edge(int u, int v, i64 cap, i64 cost) {
    g[u].push_back(es.size());
    es.emplace_back(u, v, cap, cost);
    g[v].push_back(es.size());
    es.emplace_back(v, u, 0, -cost);
  }

  bool spfa() {
    d.assign(n + 1, INF);
    in.assign(n + 1, 0);
    d[s] = 0;
    in[s] = 1;
    a[s] = INF;
    queue<int> q;
    q.push(s);
    while (!q.empty()) {
      int u = q.front();
      q.pop();
      in[u] = 0;
      for (int& i : g[u]) {
        Edge& e = es[i];
        if (e.cap && d[e.to] > d[u] + e.cost) {
          d[e.to] = d[u] + e.cost;
          p[e.to] = i;
          a[e.to] = min(a[u], e.cap);
          if (!in[e.to]) {
            q.push(e.to);
            in[e.to] = 1;
          }
        }
      }
    }
    return d[t] != INF;
  }

  void solve() {
    while (spfa()) {
      flow += a[t];
      cost += a[t] * d[t];
      int u = t;
      while (u != s) {
        es[p[u]].cap -= a[t];
        es[p[u] ^ 1].cap += a[t];
        u = es[p[u]].from;
      }
    }
  }
};
```

### Minimum Cut of Undirected Graph

```cpp
namespace stoer_wagner {
  bool vis[N], in[N];
  int g[N][N], w[N];

  void init() {
    memset(g, 0, sizeof(g));
    memset(in, 0, sizeof(in));
  }

  void add_edge(int u, int v, int w) {
    g[u][v] += w;
    g[v][u] += w;
  }

  int search(int& s, int& t) {
    memset(vis, 0, sizeof(vis));
    memset(w, 0, sizeof(w));
    int maxw, tt = n + 1;
    for (int i = 0; i < n; i++) {
      maxw = -INF;
      for (int j = 0; j < n; j++) {
        if (!in[j] && !vis[j] && w[j] > maxw) {
          maxw = w[j];
          tt = j;
        }
      }
      if (t == tt) return w[t];
      s = t; t = tt;
      vis[tt] = true;
      for (int j = 0; j < n; j++) {
        if (!in[j] && !vis[j]) {
          w[j] += g[tt][j];
        }
      }
    }
    return w[t];
  }

  int go() {
    int s, t, ans = INF;
    for (int i = 0; i < n - 1; i++) {
      s = t = -1;
      ans = min(ans, search(s, t));
      if (ans == 0) return 0;
      in[t] = true;
      for (int j = 0; j < n; j++) {
        if (!in[j]) {
          g[s][j] += g[t][j];
          g[j][s] += g[j][t];
        }
      }
    }
    return ans;
  }
}
```

### Heavy-Light Decomposition

```cpp
// weights on vertices
vector<int> g[N];
int pa[N], sz[N], dep[N], dfn[N], maxc[N], top[N], clk;

void dfs1(int u) {
  sz[u] = 1;
  maxc[u] = -1;
  int maxs = 0;
  for (int& v : g[u]) {
    if (v != pa[u]) {
      pa[v] = u;
      dep[v] = dep[u] + 1;
      dfs1(v);
      sz[u] += sz[v];
      if (umax(maxs, sz[v])) maxc[u] = v;
    }
  }
}

void dfs2(int u, int tp) {
  top[u] = tp;
  dfn[u] = ++clk;
  if (maxc[u] != -1) dfs2(maxc[u], tp);
  for (int& v : g[u]) {
    if (v != pa[u] && v != maxc[u]) {
      dfs2(v, v);
    }
  }
}

void init() {
  clk = 0;
  dep[1] = 1;
  dfs1(1);
  dfs2(1, 1);
}

i64 go(int u, int v) {
  int uu = top[u], vv = top[v];
  i64 res = 0;
  while (uu != vv) {
    if (dep[uu] < dep[vv]) {
      swap(u, v);
      swap(uu, vv);
    }
    res += segt.query(dfn[uu], dfn[u]);
    u = pa[uu];
    uu = top[u];
  }
  if (dep[u] > dep[v]) swap(u, v);
  res += segt.query(dfn[u], dfn[v]);
  return res;
}
```

### Tarjan

+ Cut Points

```cpp
int dfn[N], low[N], clk;

void init() { clk = 0; memset(dfn, 0, sizeof(dfn)); }

void tarjan(int u, int pa) {
  low[u] = dfn[u] = ++clk;
  int cc = (pa != 0);
  for (int v : g[u]) {
    if (v == pa) continue;
    if (!dfn[v]) {
      tarjan(v, u);
      low[u] = min(low[u], low[v]);
      cc += low[v] >= dfn[u];
    } else low[u] = min(low[u], dfn[v]);
  }
  if (cc > 1) // ...
}
```

+ Bridges

```cpp
int dfn[N], low[N], clk;

void init() { clk = 0; memset(dfn, 0, sizeof(dfn)); }

void tarjan(int u, int pa) {
  low[u] = dfn[u] = ++clk;
  int f = 0;
  for (int v : g[u]) {
    if (v == pa && ++f == 1) continue;
    if (!dfn[v]) {
      tarjan(v, u);
      if (low[v] > dfn[u]) // ...
      low[u] = min(low[u], low[v]);
    } else low[u] = min(low[u], dfn[v]);
  }
}
```

+ Strongly Connected Components (SCC)

```cpp
int dfn[N], low[N], clk, tot, color[N];
vector<int> scc[N];

void init() { tot = clk = 0; memset(dfn, 0, sizeof dfn); }

void tarjan(int u) {
  static int st[N], p;
  static bool in[N];
  dfn[u] = low[u] = ++clk;
  st[p++] = u;
  in[u] = true;
  for (int v : g[u]) {
    if (!dfn[v]) {
      tarjan(v);
      low[u] = min(low[u], low[v]);
    } else if (in[v]) {
      low[u] = min(low[u], dfn[v]);
    }
  }
  if (dfn[u] == low[u]) {
    ++tot;
    for (;;) {
      int x = st[--p];
      in[x] = false;
      color[x] = tot;
      scc[tot].push_back(x);
      if (x == u) break;
    }
  }
}
```

+ 2-SAT

```cpp
// N doubled
void two_sat() {
  for (int i = 1; i <= n * 2; i++) {
    if (!dfn[i]) tarjan(i);
  }
  for (int i = 1; i <= n; i++) {
    if (color[i] == color[i + n]) {
      // impossible
    }
  }
  for (int i = 1; i <= n; i++) {
    if (color[i] < color[i + n]) {
      // select
    }
  }
}
```

### Eulerian Path

+ Undirected Graph

```cpp
vector<int> euler_path(int s) {
  vector<int> path;
  stack<pair<int, int>> st;
  st.emplace(s, -1);
  while (!st.empty()) {
    auto [u, i] = st.top();
    if (g[u].empty()) {
      if (i != -1) path.push_back(i);
      st.pop();
    } else {
      i = *g[u].begin();
      int v = es[i].first ^ es[i].second ^ u;
      g[u].erase(i);
      g[v].erase(i);
      st.emplace(v, i);
    }
  }
  return path;
}
```

+ Directed Graph

```cpp
vector<int> euler_path(int s) {
  vector<int> path;
  stack<pair<int, int>> st;
  st.emplace(s, -1);
  while (!st.empty()) {
    auto [u, i] = st.top();
    if (g[u].empty()) {
      if (i != -1) path.push_back(i);
      st.pop();
    } else {
      i = *g[u].begin();
      int v = es[i].second;
      g[u].erase(i);
      st.emplace(v, i);
    }
  }
  reverse(path.begin(), path.end());
  return path;
}
```

### Dominator Tree

+ Directed Acyclic Graph (DAG)

```cpp
// rt is a point in g with in-degree 0 (may need to create a super source point)
const int LOG = 22;
int n, deg[N], dep[N], up[N][LOG];
vector<int> g[N], rg[N], dt[N];

bool topo(vector<int>& ans, int rt) {
  queue<int> q;
  q.push(rt);
  while (!q.empty()) {
    int u = q.front();
    q.pop();
    ans.push_back(u);
    for (int v : g[u]) {
      deg[v]--;
      if (deg[v] == 0) q.push(v);
    }
  }
  return ans.size() == n;
}

int lca(int u, int v) {
  if (dep[u] > dep[v]) swap(u, v);
  int t = dep[v] - dep[u];
  for (int i = 0; i < LOG; i++) {
    if ((t >> i) & 1) v = up[v][i];
  }
  if (u == v) return u;
  for (int i = LOG - 1; i >= 0; i--) {
    if (up[u][i] != up[v][i]) {
      u = up[u][i];
      v = up[v][i];
    }
  }
  return up[u][0];
}

void go(int rt) {
  vector<int> a;
  topo(a, rt);
  dep[rt] = 1;
  for (int i = 1; i < a.size(); i++) {
    int u = a[i], pa = -1;
    for (int v : rg[u]) {
      pa = (pa == -1) ? v : lca(pa, v);
    }
    dt[pa].push_back(u);
    dep[u] = dep[pa] + 1;
    up[u][0] = pa;
    for (int i = 1; i < LOG; i++) {
      up[u][i] = up[up[u][i - 1]][i - 1];
    }
  }
}
```

+ General Directed Graph

```cpp
vector<int> g[N], rg[N];
vector<int> dt[N];

namespace tl {
  int pa[N], dfn[N], clk, rdfn[N];
  int c[N], best[N], sdom[N], idom[N];

  void init(int n) {
    clk = 0;
    fill(c, c + n + 1, -1);
    fill(dfn, dfn + n + 1, 0);
    for (int i = 1; i <= n; i++) {
      dt[i].clear();
      sdom[i] = best[i] = i;
    }
  }

  void dfs(int u) {
    dfn[u] = ++clk;
    rdfn[clk] = u;
    for (int& v : g[u]) {
      if (!dfn[v]) {
        pa[v] = u;
        dfs(v);
      }
    }
  }

  int fix(int x) {
    if (c[x] == -1) return x;
    int& f = c[x], rt = fix(f);
    if (dfn[sdom[best[x]]] > dfn[sdom[best[f]]]) best[x] = best[f];
    return f = rt;
  }

  void go(int rt) {
    dfs(rt);
    for (int i = clk; i > 1; i--) {
      int x = rdfn[i], mn = clk + 1;
      for (int& u : rg[x]) {
        if (!dfn[u]) continue; // may not reach all vertices
        fix(u);
        mn = min(mn, dfn[sdom[best[u]]]);
      }
      c[x] = pa[x];
      dt[sdom[x] = rdfn[mn]].push_back(x);
      x = rdfn[i - 1];
      for (int& u: dt[x]) {
        fix(u);
        idom[u] = (sdom[best[u]] == x) ? x : best[u];
      }
      dt[x].clear();
    }
    for (int i = 2; i <= clk; i++) {
      int u = rdfn[i];
      if (idom[u] != sdom[u]) idom[u] = idom[idom[u]];
      dt[idom[u]].push_back(u);
    }
  }
}
```
