## 字符串

### 哈希

```cpp
// open hack
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());

i64 random_prime(i64 l, i64 r) {
  for (;;) {
    i64 x = rng() % (r - l + 1) + l;
    if (isprime(x)) return x;
  }
}

const int x = 135, p = 1e9 + 9;

i64 xp[N];

void init_xp() {
  xp[0] = 1;
  for (int i = 1; i < N; i++) {
    xp[i] = xp[i - 1] * x % p;
  }
}

struct Hash {
  vector<i64> h;

  Hash() : h(1) {}

  void add(const string &s) {
    i64 res = h.back();
    for (char c : s) {
      res = (res * x + c) % p;
      h.push_back(res);
    }
  }

  // 0-indexed, [l, r]
  i64 get(int l, int r) {
    r++;
    return (h[r] - h[l] * xp[r - l] % p + p) % p;
  }
};
```

+ 双哈希

```cpp
const int x = 135, p1 = 1e9 + 7, p2 = 1e9 + 9;
const ull mask32 = ~(0u);

ull xp1[N], xp2[N];

void init_xp() {
  xp1[0] = xp2[0] = 1;
  for (int i = 1; i < N; i++) {
    xp1[i] = xp1[i - 1] * x % p1;
    xp2[i] = xp2[i - 1] * x % p2;
  }
}

struct Hash {
  vector<ull> h;

  Hash() : h(1) {}

  void add(const string& s) {
    ull res1 = h.back() >> 32;
    ull res2 = h.back() & mask32;
    for (char c : s) {
      res1 = (res1 * x + c) % p1;
      res2 = (res2 * x + c) % p2;
      h.push_back((res1 << 32) | res2);
    }
  }

  // 0-indexed, [l, r]
  ull get(int l, int r) {
    r++;
    int len = r - l;
    ull l1 = h[l] >> 32, r1 = h[r] >> 32;
    ull l2 = h[l] & mask32, r2 = h[r] & mask32;
    ull res1 = (r1 - l1 * xp1[len] % p1 + p1) % p1;
    ull res2 = (r2 - l2 * xp2[len] % p2 + p2) % p2;
    return (res1 << 32) | res2;
  }
};
```

+ 二维哈希

```cpp
const i64 basex = 239, basey = 241, p = 998244353;

i64 pwx[N], pwy[N];

void init_xp() {
  pwx[0] = pwy[0] = 1;
  for (int i = 1; i < N; i++) {
    pwx[i] = pwx[i - 1] * basex % p;
    pwy[i] = pwy[i - 1] * basey % p;
  }
}

struct Hash2D {
  vector<vector<i64> > h;

  Hash2D(const vector<vector<int> >& a, int n, int m) : h(n + 1, vector<i64>(m + 1)) {
    for (int i = 0; i < n; i++) {
      i64 s = 0;
      for (int j = 0; j < m; j++) {
        s = (s * basey + a[i][j] + 1) % p;
        h[i + 1][j + 1] = (h[i][j + 1] * basex + s) % p;
      }
    }
  }

  i64 get(int x, int y, int xx, int yy) {
    ++xx; ++yy;
    int dx = xx - x, dy = yy - y;
    i64 res = h[xx][yy]
      - h[x][yy] * pwx[dx]
      - h[xx][y] * pwy[dy]
      + h[x][y] * pwx[dx] % p * pwy[dy];
    return (res % p + p) % p;
  }
};
```

+ 动态哈希

```cpp
const i64 MAGIC = 479;
const i64 IMAGIC = 628392490;
const i64 P = 1e9 + 9;
const int N = 1e5 + 5;

i64 pw[N], inv[N];

void init() {
  pw[0] = inv[0] = 1;
  for (int i = 1; i < N; i++) {
    pw[i] = (pw[i - 1] * MAGIC) % P;
    inv[i] = (inv[i - 1] * IMAGIC) % P;
  }
}

// 1-indexed, [l, r]
struct fenwick {
  int n;
  vector<i64> t;
  fenwick(int n) : n(n), t(n + 1) {}
  void add(int p, i64 x) {
    for (; p <= n; p += p & -p) (t[p] += x) %= P;
  }
  i64 get(int p) {
    i64 a = 0;
    for (; p > 0; p -= p & -p) (a += t[p]) %= P;
    return a;
  }
  void set(int p, i64 x) { add(p, (x - query(p, p) + P) * pw[p] % P); }
  i64 query(int l, int r) { return (get(r) - get(l - 1) + P) * inv[l] % P; }
};
```

### Manacher

```cpp
// "aba" => "#a#b#a#"
struct Manacher {
  vector<int> d;

  Manacher(const string& s) {
    string t = "#";
    for (int i = 0; i < s.size(); i++) {
      t.push_back(s[i]);
      t.push_back('#');
    }
    int n = t.size();
    d.resize(n);
    for (int i = 0, l = 0, r = -1; i < n; i++) {
      int k = (i > r) ? 1 : min(d[l + r - i], r - i);
      while (i - k >= 0 && i + k < n && t[i - k] == t[i + k]) k++;
      d[i] = --k;
      if (i + k > r) {
        l = i - k;
        r = i + k;
      }
    }
  }

  // 0-indexed [l, r]
  bool is_p(int l, int r) { return d[l + r + 1] >= r - l + 1; }
};
```

### KMP

```cpp
// 前缀函数（每一个前缀的最长公共前后缀）
vector<int> get_pi(const string& s) {
  int n = s.size();
  vector<int> a(n);
  for (int i = 1, j = 0; i < n; i++) {
    while (j && s[j] != s[i]) j = a[j - 1];
    if (s[j] == s[i]) j++;
    a[i] = j;
  }
  return a;
}

struct KMP {
  string s;
  vector<int> a;

  KMP(const string& s) : s(s), a(get_pi(s)) {}

  vector<int> find_in(const string& t) {
    vector<int> res;
    for (int i = 0, j = 0; i < t.size(); i++) {
      while (j && s[j] != t[i]) j = a[j - 1];
      if (s[j] == t[i]) j++;
      if (j == s.size()) {
        res.push_back(i - j + 1);
        j = a[j - 1]; // 允许重叠匹配 j = 0 不允许
      }
    }
    return res;
  }
};

// Z函数（每一个后缀和该字符串的最长公共前缀）
vector<int> get_z(const string& s) {
  int n = s.size(), l = 0, r = 0;
  vector<int> z(n);
  for (int i = 1; i < n; i++) {
    if (i <= r) z[i] = min(r - i + 1, z[i - l]);
    while (i + z[i] < n && s[z[i]] == s[i + z[i]]) z[i]++;
    if (i + z[i] - 1 > r) {
      l = i;
      r = i + z[i] - 1;
    }
  }
  return z;
}
```

### Lyndon 分解

```cpp
vector<string> duval(const string& s) {
  int n = s.size(), i = 0;
  vector<string> d;
  while (i < n) {
    int j = i + 1, k = i;
    while (j < n && s[k] <= s[j]) {
      if (s[k] < s[j]) k = i;
      else k++;
      j++;
    }
    while (i <= k) {
      d.push_back(s.substr(i, j - k));
      i += j - k;
    }
  }
  return d;
}
```

### 最小表示法

```cpp
int get(const string& s) {
  int k = 0, i = 0, j = 1, n = s.size();
  while (k < n && i < n && j < n) {
    if (s[(i + k) % n] == s[(j + k) % n]) {
      k++;
    } else {
      s[(i + k) % n] > s[(j + k) % n] ? i = i + k + 1 : j = j + k + 1;
      if (i == j) i++;
      k = 0;
    }
  }
  return min(i, j);
}
```

### Trie

```cpp
// 01 Trie
struct Trie {
  int t[31 * N][2], sz;

  void init() {
    memset(t, 0, 2 * (sz + 2) * sizeof(int));
    sz = 1;
  }

  void insert(int x) {
    int p = 0;
    for (int i = 30; i >= 0; i--) {
      bool d = (x >> i) & 1;
      if (!t[p][d]) t[p][d] = sz++;
      p = t[p][d];
    }
  }
};

// 正常Trie
struct Trie {
  int t[N][26], sz, cnt[N];

  void init() {
    memset(t, 0, 26 * (sz + 2) * sizeof(int));
    memset(cnt, 0, (sz + 2) * sizeof(int));
    sz = 1;
  }

  void insert(const string& s) {
    int p = 0;
    for (char c : s) {
      int d = c - 'a';
      if (!t[p][d]) t[p][d] = sz++;
      p = t[p][d];
    }
    cnt[p]++;
  }
};
```

### AC 自动机

```cpp
struct ACA {
  int t[N][26], sz, fail[N], nxt[N], cnt[N];

  void init() {
    memset(t, 0, 26 * (sz + 2) * sizeof(int));
    memset(fail, 0, (sz + 2) * sizeof(int));
    memset(nxt, 0, (sz + 2) * sizeof(int));
    memset(cnt, 0, (sz + 2) * sizeof(int));
    sz = 1;
  }

  void insert(const string& s) {
    int p = 0;
    for (char c : s) {
      int d = c - 'a';
      if (!t[p][d]) t[p][d] = sz++;
      p = t[p][d];
    }
    cnt[p]++;
  }

  void build() {
    queue<int> q;
    for (int i = 0; i < 26; i++) {
      if (t[0][i]) q.push(t[0][i]);
    }
    while (!q.empty()) {
      int u = q.front();
      q.pop();
      for (int i = 0; i < 26; i++) {
        int& v = t[u][i];
        if (v) {
          fail[v] = t[fail[u]][i];
          nxt[v] = cnt[fail[v]] ? fail[v] : nxt[fail[v]];
          q.push(v);
        } else {
          v = t[fail[u]][i];
        }
      }
    }
  }
};
```

### 回文自动机

```cpp
// WindJ0Y
struct Palindromic_Tree {
  static constexpr int N = 300005;

  int next[N][26]; // next指针，next指针和字典树类似，指向的串为当前串两端加上同一个字符构成
  int fail[N]; // fail指针，失配后跳转到fail指针指向的节点
  int cnt[N]; // 表示节点i表示的本质不同的串的个数 after count()
  int num[N]; // 表示以节点i表示的最长回文串的最右端点为回文串结尾的回文串个数。
  int len[N]; // len[i]表示节点i表示的回文串的长度
  int lcnt[N];
  int S[N]; // 存放添加的字符
  int last; // 指向上一个字符所在的节点，方便下一次add
  int n; // 字符数组指针
  int p; // 节点指针

  int newnode(int l, int vc) { // 新建节点
    for (int i = 0; i < 26; ++i) next[p][i] = 0;
    cnt[p] = 0;
    num[p] = 0;
    len[p] = l;
    lcnt[p] = vc;
    return p++;
  }

  void init() { // 初始化
    p = 0;
    newnode(0, 0);
    newnode(-1, 0);
    last = 0;
    n = 0;
    S[n] = -1; // 开头放一个字符集中没有的字符，减少特判
    fail[0] = 1;
  }

  int get_fail(int x) { // 和KMP一样，失配后找一个尽量最长的
    while (S[n - len[x] - 1] != S[n]) x = fail[x];
    return x;
  }

  void add(int c) {
    S[++n] = c;
    int cur = get_fail(last); // 通过上一个回文串找这个回文串的匹配位置
    if (!next[cur][c]) { // 如果这个回文串没有出现过，说明出现了一个新的本质不同的回文串
      int now = newnode(len[cur] + 2, lcnt[cur] | (1 << c)); // 新建节点
      fail[now] = next[get_fail(fail[cur])][c]; // 和AC自动机一样建立fail指针，以便失配后跳转
      next[cur][c] = now;
      num[now] = num[fail[now]] + 1;
    }
    last = next[cur][c];
    cnt[last]++;
  }

  void count() {
    for (int i = p - 1; i >= 0; --i) cnt[fail[i]] += cnt[i];
    // 父亲累加儿子的cnt，因为如果fail[v]=u，则u一定是v的子回文串
  }
} pt;
```

### 后缀自动机

```cpp
// 下标从 1 开始
// rsort 中的数组 a 是拓扑序 [1, sz)
struct SAM {
  static constexpr int M = N << 1;
  int t[M][26], len[M], fa[M], sz = 2, last = 1;
  void init() {
    memset(t, 0, (sz + 2) * sizeof t[0]);
    sz = 2;
    last = 1;
  }
  void ins(int ch) {
    int p = last, np = last = sz++;
    len[np] = len[p] + 1;
    for (; p && !t[p][ch]; p = fa[p]) t[p][ch] = np;
    if (!p) {
      fa[np] = 1;
      return;
    }
    int q = t[p][ch];
    if (len[q] == len[p] + 1) {
      fa[np] = q;
    } else {
      int nq = sz++;
      len[nq] = len[p] + 1;
      memcpy(t[nq], t[q], sizeof t[0]);
      fa[nq] = fa[q];
      fa[np] = fa[q] = nq;
      for (; t[p][ch] == q; p = fa[p]) t[p][ch] = nq;
    }
  }

  int c[M] = {1}, a[M];
  void rsort() {
    for (int i = 1; i < sz; ++i) c[i] = 0;
    for (int i = 1; i < sz; ++i) c[len[i]]++;
    for (int i = 1; i < sz; ++i) c[i] += c[i - 1];
    for (int i = 1; i < sz; ++i) a[--c[len[i]]] = i;
  }
};
```

+ 广义后缀自动机（在线版）

```cpp
// 插入新串前 置 last 为 1
struct SAM {
  static constexpr int M = N << 1;
  int t[M][26], len[M], fa[M], sz = 2, last = 1;
  void init() {
    memset(t, 0, (sz + 2) * sizeof t[0]);
    sz = 2;
    last = 1;
  }
  void ins(int ch) {
    int p = last, np = 0, nq = 0, q = -1;
    if (!t[p][ch]) {
      np = sz++;
      len[np] = len[p] + 1;
      for (; p && !t[p][ch]; p = fa[p]) t[p][ch] = np;
    }
    if (!p) {
      fa[np] = 1;
    } else {
      q = t[p][ch];
      if (len[p] + 1 == len[q]) {
        fa[np] = q;
      } else {
        nq = sz++;
        len[nq] = len[p] + 1;
        memcpy(t[nq], t[q], sizeof t[0]);
        fa[nq] = fa[q];
        fa[np] = fa[q] = nq;
        for (; t[p][ch] == q; p = fa[p]) t[p][ch] = nq;
      }
    }
    last = np ? np : nq ? nq : q;
  }
};
```

### 后缀数组

```cpp
// 下标从1开始
// sa[i]: 排名为i的后缀位置
// rk[i]: 第i个后缀的排名
// ht[i]: LCP(sa[i], sa[i - 1])
struct SA {
  int n, m;
  vector<int> a, d, sa, rk, ht;

  void rsort() {
    vector<int> c(m + 1);
    for (int i = 1; i <= n; i++) c[rk[d[i]]]++;
    for (int i = 1; i <= m; i++) c[i] += c[i - 1];
    for (int i = n; i; i--) sa[c[rk[d[i]]]--] = d[i];
  }

  SA(const string& s) : n(s.size()), m(128), a(n + 1), d(n + 1), sa(n + 1), rk(n + 1), ht(n + 1) {
    for (int i = 1; i <= n; i++) { rk[i] = a[i] = s[i - 1]; d[i] = i; }
    rsort();
    for (int j = 1, i, k; k < n; m = k, j <<= 1) {
      for (i = n - j + 1, k = 0; i <= n; i++) d[++k] = i;
      for (i = 1; i <= n; i++) if (sa[i] > j) d[++k] = sa[i] - j;
      rsort(); swap(rk, d); rk[sa[1]] = k = 1;
      for (i = 2; i <= n; i++) {
        rk[sa[i]] = (d[sa[i]] == d[sa[i - 1]] && d[sa[i] + j] == d[sa[i - 1] + j]) ? k : ++k;
      }
    }
    int j, k = 0;
    for (int i = 1; i <= n; ht[rk[i++]] = k) {
      for (k ? k-- : k, j = sa[rk[i] - 1]; a[i + k] == a[j + k]; ++k);
    }
  }
};
```
