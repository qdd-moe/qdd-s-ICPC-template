## Strings

### Hash

```cpp
// open hack
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());

struct Hash {
  static const i64 md = (1LL << 61) - 1;
  static i64 step;
  static vector<i64> pw;

  static i64 mul(i64 a, i64 b) { return (i64)(__int128(a) * b % md); }

  static void init(int N) {
    pw.resize(N + 1);
    pw[0] = 1;
    for (int i = 1; i <= N; i++) pw[i] = mul(pw[i - 1], step);
  }

  vector<i64> h;

  template <class T>
  Hash(const T& s) {
    int n = s.size();
    h.resize(n + 1);
    for (int i = 0; i < n; i++) h[i + 1] = (mul(h[i], step) + s[i]) % md;
  }

  i64 operator()(int l, int r) { return (h[r + 1] - mul(h[l], pw[r - l + 1]) + md) % md; }
};

i64 Hash::step = uniform_int_distribution<i64>(256, md - 1)(rng);
vector<i64> Hash::pw;
int _ = (Hash::init(1e6), 0);
```

+ 2D Hash

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

+ Dynamic Hash

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

+ SplitMix

```cpp
// tourist
mt19937_64 rng((unsigned int)chrono::steady_clock::now().time_since_epoch().count());

struct hash61 {
  static const uint64_t md = (1LL << 61) - 1;
  static uint64_t step;
  static vector<uint64_t> pw;

  uint64_t addmod(uint64_t a, uint64_t b) const {
    a += b;
    if (a >= md) a -= md;
    return a;
  }

  uint64_t submod(uint64_t a, uint64_t b) const {
    a += md - b;
    if (a >= md) a -= md;
    return a;
  }

  uint64_t mulmod(uint64_t a, uint64_t b) const {
    uint64_t l1 = (uint32_t)a, h1 = a >> 32, l2 = (uint32_t)b, h2 = b >> 32;
    uint64_t l = l1 * l2, m = l1 * h2 + l2 * h1, h = h1 * h2;
    uint64_t ret = (l & md) + (l >> 61) + (h << 3) + (m >> 29) + (m << 35 >> 3) + 1;
    ret = (ret & md) + (ret >> 61);
    ret = (ret & md) + (ret >> 61);
    return ret - 1;
  }

  void ensure_pw(int sz) {
    int cur = (int)pw.size();
    if (cur < sz) {
      pw.resize(sz);
      for (int i = cur; i < sz; i++) {
        pw[i] = mulmod(pw[i - 1], step);
      }
    }
  }

  vector<uint64_t> pref;
  int n;

  template <typename T>
  hash61(const T& s) {
    n = (int)s.size();
    ensure_pw(n + 1);
    pref.resize(n + 1);
    pref[0] = 1;
    for (int i = 0; i < n; i++) {
      pref[i + 1] = addmod(mulmod(pref[i], step), s[i]);
    }
  }

  inline uint64_t operator()(const int from, const int to) const {
    assert(0 <= from && from <= to && to <= n - 1);
    return submod(pref[to + 1], mulmod(pref[from], pw[to - from + 1]));
  }
};

uint64_t hash61::step = (md >> 2) + rng() % (md >> 1);
vector<uint64_t> hash61::pw = vector<uint64_t>(1, 1);
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
// prefix function (the longest common prefix and suffix for each prefix)
// the smallest cycle of [0...i] is i + 1 - a[i]
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
        j = a[j - 1]; // Allowing overlapping matches j = 0 not allowed
      }
    }
    return res;
  }
};

// Z function, z[i] = LCP(s, s[i:])
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

### Lyndon Decomposition

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

### Lexicographically Minimal String Rotation

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

// Normal Trie
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

### Aho-Corasick Automaton

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

### Palindromic Tree

```cpp
// WindJ0Y
struct Palindromic_Tree {
  static constexpr int N = 300005;

  int next[N][26]; // next pointer, similar to trie, points to the string formed by adding the same character at the beginning and end of the current string
  int fail[N]; // fail pointer, jumps to the node pointed to by fail pointer after mismatch
  int cnt[N]; // represents the number of different essential strings represented by node i after count()
  int num[N]; // represents the number of palindrome strings with the last character of the palindrome string represented by node i as the end of the palindrome string.
  int len[N]; // len[i] represents the length of the palindrome string represented by node i
  int lcnt[N];
  int S[N]; // Store the added characters
  int last; // Points to the node where the last character is located for easy addition next time
  int n; // Character array pointer
  int p; // Node pointer

  int newnode(int l, int vc) { // Create a new node
    for (int i = 0; i < 26; ++i) next[p][i] = 0;
    cnt[p] = 0;
    num[p] = 0;
    len[p] = l;
    lcnt[p] = vc;
    return p++;
  }

  void init() { // Initialize
    p = 0;
    newnode(0, 0);
    newnode(-1, 0);
    last = 0;
    n = 0;
    S[n] = -1; // Put a character not in the character set at the beginning to reduce special cases
    fail[0] = 1;
  }

  int get_fail(int x) { // Similar to KMP, find the longest one after mismatch
    while (S[n - len[x] - 1] != S[n]) x = fail[x];
    return x;
  }

  void add(int c) {
    S[++n] = c;
    int cur = get_fail(last); // Find the matching position of this palindrome string through the last palindrome string
    if (!next[cur][c]) { // If this palindrome string has not appeared, it means a new essential different palindrome string has appeared
      int now = newnode(len[cur] + 2, lcnt[cur] | (1 << c)); // Create a new node
      fail[now] = next[get_fail(fail[cur])][c]; // Establish a fail pointer as in AC automation, so as to jump after mismatch
      next[cur][c] = now;
      num[now] = num[fail[now]] + 1;
    }
    last = next[cur][c];
    cnt[last]++;
  }

  void count() {
    for (int i = p - 1; i >= 0; --i) cnt[fail[i]] += cnt[i];
    // Parent node accumulates child node's cnt, because if fail[v]=u, then u must be v's child palindrome string
  }
} pt;
```

### Suffix Automaton

```cpp
// 1-indexed
// the array a in rsort is the topological order [1, sz)
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

+ Multiple Strings, Online Construction

```cpp
// set last to 1 before inserting a new string
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

### Suffix Array

```cpp
// 0-indexed
// sa[i]: position of the suffix with rank i
// rk[i]: rank of the i-th suffix
// lc[i]: LCP(sa[i], sa[i + 1])
struct SuffixArray {
  int n;
  vector<int> sa, rk, lc;
  SuffixArray(const string& s) {
    n = s.length();
    sa.resize(n);
    lc.resize(n - 1);
    rk.resize(n);
    iota(sa.begin(), sa.end(), 0);
    sort(sa.begin(), sa.end(), [&](int a, int b) { return s[a] < s[b]; });
    rk[sa[0]] = 0;
    for (int i = 1; i < n; ++i) rk[sa[i]] = rk[sa[i - 1]] + (s[sa[i]] != s[sa[i - 1]]);
    int k = 1;
    vector<int> tmp, cnt(n);
    tmp.reserve(n);
    while (rk[sa[n - 1]] < n - 1) {
      tmp.clear();
      for (int i = 0; i < k; ++i) tmp.push_back(n - k + i);
      for (auto i : sa)
        if (i >= k) tmp.push_back(i - k);
      fill(cnt.begin(), cnt.end(), 0);
      for (int i = 0; i < n; ++i) ++cnt[rk[i]];
      for (int i = 1; i < n; ++i) cnt[i] += cnt[i - 1];
      for (int i = n - 1; i >= 0; --i) sa[--cnt[rk[tmp[i]]]] = tmp[i];
      swap(rk, tmp);
      rk[sa[0]] = 0;
      for (int i = 1; i < n; ++i)
        rk[sa[i]] = rk[sa[i - 1]] + (tmp[sa[i - 1]] < tmp[sa[i]] || sa[i - 1] + k == n || tmp[sa[i - 1] + k] < tmp[sa[i] + k]);
      k *= 2;
    }
    for (int i = 0, j = 0; i < n; ++i) {
      if (rk[i] == 0) {
        j = 0;
      } else {
        for (j -= j > 0; i + j < n && sa[rk[i] - 1] + j < n && s[i + j] == s[sa[rk[i] - 1] + j];) ++j;
        lc[rk[i] - 1] = j;
      }
    }
  }
};
```
