## Math

### GCD & LCM

```cpp
i64 gcd(i64 a, i64 b) { return b ? gcd(b, a % b) : a; }
i64 lcm(i64 a, i64 b) { return a / gcd(a, b) * b; }
```

### Safe Multiplication & Fast Exponentiation

```cpp
// Use if modulus exceeds INT_MAX
i64 mul(i64 a, i64 b, i64 p) {
  i64 ans = 0;
  for (a %= p; b; b >>= 1, a = (a << 1) % p)
    if (b & 1) ans = (ans + a) % p;
  return ans;
}

// O(1)
i64 mul(i64 a, i64 b, i64 p) {
  return (i64)(__int128(a) * b % p);
}

i64 qk(i64 a, i64 b, i64 p) {
  i64 ans = 1 % p;
  for (a %= p; b; b >>= 1, a = a * a % p)
    if (b & 1) ans = ans * a % p;
  return ans;
}

// If modulus exceeds INT_MAX
i64 qk(i64 a, i64 b, i64 p) {
  i64 ans = 1 % p;
  for (a %= p; b; b >>= 1, a = mul(a, a, p))
    if (b & 1) ans = mul(ans, a, p);
  return ans;
}

// Decimal fast exponentiation
i64 qk(i64 a, const string& b, i64 p) {
  i64 ans = 1;
  for (int i = b.size() - 1; i >= 0; i--) {
    ans = ans * qk(a, b[i] - '0', p) % p;
    a = qk(a, 10, p);
  }
  return ans;
}
```

### Matrix Fast Exponentiation

```cpp
const int M_SZ = 3;

using Mat = array<array<i64, M_SZ>, M_SZ>;
using Vec = array<i64, M_SZ>;

#define rep2 for (int i = 0; i < M_SZ; i++) for (int j = 0; j < M_SZ; j++)

void zero(Mat& a) { rep2 a[i][j] = 0; }
void one(Mat& a) { rep2 a[i][j] = (i == j); }

Vec mul(const Mat& a, const Vec& b, i64 p) {
  Vec ans;
  fill(ans.begin(), ans.end(), 0);
  rep2 { (ans[i] += a[i][j] * b[j]) %= p; }
  return ans;
}

Mat mul(const Mat& a, const Mat& b, i64 p) {
  Mat ans; zero(ans);
  rep2 if (a[i][j]) for (int k = 0; k < M_SZ; k++) {
    (ans[i][k] += a[i][j] * b[j][k]) %= p;
  }
  return ans;
}

Mat qk(Mat a, i64 b, i64 p) {
  Mat ans; one(ans);
  for (; b; b >>= 1) {
    if (b & 1) ans = mul(a, ans, p);
    a = mul(a, a, p);
  }
  return ans;
}

// Decimal fast exponentiation
Mat qk(Mat a, const string& b, i64 p) {
  Mat ans; one(ans);
  for (int i = b.size() - 1; i >= 0; i--) {
    ans = mul(qk(a, b[i] - '0', p), ans, p);
    a = qk(a, 10, p);
  }
  return ans;
}

#undef rep2
```

### Prime Test

```cpp
bool isprime(int x) {
  if (x < 2) return false;
  for (int i = 2; i * i <= x; i++) {
    if (x % i == 0) return false;
  }
  return true;
}
```

### Sieves

```cpp
struct sieve {
  vector<bool> vis;
  vector<int> prime, spf, phi, mu;

  void init(int n) {
    vis.assign(n + 1, false);
    spf.resize(n + 1);
    phi.resize(n + 1);
    mu.resize(n + 1);
    spf[1] = phi[1] = mu[1] = 1;
    for (int i = 2; i <= n; i++) {
      if (!vis[i]) {
        prime.push_back(i);
        spf[i] = i;
        phi[i] = i - 1;
        mu[i] = -1;
      }
      for (int p : prime) {
        int d = i * p;
        if (d > n) break;
        vis[d] = true;
        spf[d] = p;
        if (i % p == 0) {
          phi[d] = phi[i] * p;
          mu[d] = 0;
          break;
        } else {
          phi[d] = phi[i] * (p - 1);
          mu[d] = -mu[i];
        }
      }
    }
  }
};
```

### Interval Sieve

```cpp
// a, b <= 1e13, b - a <= 1e6
bool vis_small[N], vis_big[N];
i64 prime[N];
int tot = 0;

void get_prime(i64 a, i64 b) {
  i64 c = ceil(sqrt(b));
  for (i64 i = 2; i <= c; i++) {
    if (!vis_small[i]) {
      for (i64 j = i * i; j <= c; j += i) {
        vis_small[j] = 1;
      }
      for (i64 j = max(i, (a + i - 1) / i) * i; j <= b; j += i) {
        vis_big[j - a] = 1;
      }
    }
  }
  for (int i = max(0LL, 2 - a); i <= b - a; i++) {
    if (!vis_big[i]) prime[tot++] = i + a;
  }
}
```

### Find Factors

```cpp
// O(sqrt(n))
vector<int> getf(int x) {
  vector<int> v;
  for (int i = 1; i * i <= x; i++) {
    if (x % i == 0) {
      v.push_back(i);
      if (x / i != i) v.push_back(x / i);
    }
  }
  sort(v.begin(), v.end());
  return v;
}
```

### Prime Factorization

```cpp
// O(sqrt(n))
vector<pair<int, int>> getf(int x) {
  vector<pair<int, int>> v;
  for (int i = 2; i * i <= x; i++) {
    if (x % i == 0) {
      int cnt = 0;
      while (x % i == 0) x /= i, cnt++;
      v.emplace_back(i, cnt);
    }
  }
  if (x != 1) v.emplace_back(x, 1);
  return v;
}

// Prerequisite: Sieves
// O(logn)
vector<pair<int, int>> getf(int x) {
  vector<pair<int, int>> v;
  while (x > 1) {
    int p = spf[x];
    int cnt = 0;
    while (x % p == 0) x /= p, cnt++;
    v.emplace_back(p, cnt);
  }
  return v;
}
```

### Miller & Pollard

```cpp
i64 mul(i64 a, i64 b, i64 p) { return (i64)(__int128(a) * b % p); }

i64 qk(i64 a, i64 b, i64 p) {
  i64 ans = 1 % p;
  for (a %= p; b; b >>= 1, a = mul(a, a, p))
    if (b & 1) ans = mul(ans, a, p);
  return ans;
}

// O(logn)
// Only need to check 2, 7, 61 for 32-bit ints
bool isprime(i64 n) {
  if (n < 3) return n == 2;
  if (!(n & 1)) return false;
  i64 d = n - 1, r = 0;
  while (!(d & 1)) d >>= 1, r++;
  static vector<i64> A = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
  for (i64 a : A) {
    i64 t = qk(a, d, n);
    if (t <= 1 || t == n - 1) continue;
    for (int i = 0; i < r; i++) {
      t = mul(t, t, n);
      if (t == 1) return false;
      if (t == n - 1) break;
    }
    if (t != 1 && t != n - 1) return false;
  }
  return true;
}

mt19937_64 rng(42);

i64 pollard_rho(i64 n, i64 c) {
  i64 x = rng() % (n - 1) + 1, y = x;
  auto f = [&](i64 v) {
    i64 t = mul(v, v, n) + c;
    return t < n ? t : t - n;
  };
  for (;;) {
    x = f(x); y = f(f(y));
    if (x == y) return n;
    i64 d = gcd(abs(x - y), n);
    if (d != 1) return d;
  }
}

vector<pair<i64, int>> getf(i64 x) {
  if (x <= 1) return {};
  map<i64, int> mp;
  function<void(i64)> f = [&](i64 n) {
    if (n == 4) { mp[2] += 2; return; }
    if (isprime(n)) { mp[n]++; return; }
    i64 p = n, c = 19260817;
    while (p == n) p = pollard_rho(n, --c);
    f(p); f(n / p);
  };
  f(x);
  return vector<pair<i64, int>>(mp.begin(), mp.end());
}
```

### Euler's Totient Function

```cpp
// Prerequisite: Prime Factorization
int phi(int x) {
  int ret = x;
  for (auto& [f, _] : getf(x)) ret = ret / f * (f - 1);
  return ret;
}
```

### Extended Euclidean Algorithm

```cpp
// ax + by = gcd(a, b)
// return {gcd, x, y}
array<i64, 3> exgcd(i64 a, i64 b) {
  if (b == 0) return {a, 1, 0};
  auto [d, x, y] = exgcd(b, a % b);
  return {d, y, x - a / b * y};
}
```

### Euclidean-like Functions

```cpp
// f(a,b,c,n) = ∑(i=[0,n]) (ai+b)/c
// g(a,b,c,n) = ∑(i=[0,n]) i*((ai+b)/c)
// h(a,b,c,n) = ∑(i=[0,n]) ((ai+b)/c)^2
i64 f(i64 a, i64 b, i64 c, i64 n);
i64 g(i64 a, i64 b, i64 c, i64 n);
i64 h(i64 a, i64 b, i64 c, i64 n);

i64 f(i64 a, i64 b, i64 c, i64 n) {
  if (n < 0) return 0;
  i64 m = (a * n + b) / c;
  if (a >= c || b >= c) {
    return (a / c) * n * (n + 1) / 2
    + (b / c) * (n + 1)
    + f(a % c, b % c, c, n);
  } else {
    return n * m - f(c, c - b - 1, a, m - 1);
  }
}

i64 g(i64 a, i64 b, i64 c, i64 n) {
  if (n < 0) return 0;
  i64 m = (a * n + b) / c;
  if (a >= c || b >= c) {
    return (a / c) * n * (n + 1) * (2 * n + 1) / 6
    + (b / c) * n * (n + 1) / 2
    + g(a % c, b % c, c, n);
  } else {
    return (n * (n + 1) * m
    - f(c, c - b - 1, a, m - 1)
    - h(c, c - b - 1, a, m - 1)) / 2;
  }
}

i64 h(i64 a, i64 b, i64 c, i64 n) {
  if (n < 0) return 0;
  i64 m = (a * n + b) / c;
  if (a >= c || b >= c) {
    return (a / c) * (a / c) * n * (n + 1) * (2 * n + 1) / 6
    + (b / c) * (b / c) * (n + 1)
    + (a / c) * (b / c) * n * (n + 1)
    + h(a % c, b % c, c, n)
    + 2 * (a / c) * g(a % c, b % c, c, n)
    + 2 * (b / c) * f(a % c, b % c, c, n);
  } else {
    return n * m * (m + 1)
    - 2 * g(c, c - b - 1, a, m - 1)
    - 2 * f(c, c - b - 1, a, m - 1)
    - f(a, b, c, n);
  }
}
```

### Modular Inverse

```cpp
i64 inv(i64 x) { return qk(x, P - 2, P); }

// EXGCD
// Inverse exists iff gcd(a, p) = 1
i64 inv(i64 a, i64 p) {
  auto [d, x, y] = exgcd(a, p);
  if (d == 1) return (x % p + p) % p;
  return -1;
}

// Inverse Table
i64 inv[N];

void init_inv() {
  inv[1] = 1;
  for (int i = 2; i < N; i++) {
    inv[i] = 1LL * (P - P / i) * inv[P % i] % P;
  }
}
```

### Binomial Coefficient

```cpp
// Pascal's Triangle
i64 C[N][N];

void initC() {
  C[0][0] = 1;
  for (int i = 1; i < N; i++) {
    C[i][0] = 1;
    for (int j = 1; j <= i; j++) {
      C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]) % P;
    }
  }
}

// Fast Binomial Coefficient
// For combination with repetition, N needs to be 2x
i64 fac[N], ifac[N];

void init_inv() {
  fac[0] = 1;
  for (int i = 1; i < N; i++) {
    fac[i] = fac[i - 1] * i % P;
  }
  ifac[N - 1] = qk(fac[N - 1], P - 2, P);
  for (int i = N - 2; i >= 0; i--) {
    ifac[i] = ifac[i + 1] * (i + 1) % P;
  }
}

i64 C(int n, int m) {
  if (n < m || m < 0) return 0;
  return fac[n] * ifac[m] % P * ifac[n - m] % P;
}

// Lucas Theorem
i64 C(i64 n, i64 m) {
  if (n < m || m < 0) return 0;
  if (n < P && m < P) return fac[n] * ifac[m] % P * ifac[n - m] % P;
  return C(n / P, m / P) * C(n % P, m % P) % P;
}

// Combination with repetition
i64 H(int n, int m) { return C(n + m - 1, m); }
```

### Cantor Expansion

```cpp
// Requires pre-processing of factorials
int cantor(vector<int>& s) {
  int n = s.size(), ans = 0;
  for (int i = 0; i < n - 1; i++) {
    int cnt = 0;
    for (int j = i + 1; j < n; j++) {
      if (s[j] < s[i]) cnt++;
    }
    ans += cnt * fac[n - i - 1];
  }
  return ans + 1;
}

vector<int> inv_cantor(int x, int n) {
  x--;
  vector<int> ans(n), rk(n);
  iota(rk.begin(), rk.end(), 1);
  for (int i = 0; i < n; i++) {
    int t = x / fac[n - i - 1];
    x %= fac[n - i - 1];
    ans[i] = rk[t];
    for (int j = t; rk[j] < n; j++) {
      rk[j] = rk[j + 1];
    }
  }
  return ans;
}
```

### Gauss-Jordan elimination

+ Floating-point version

```cpp
// n: number of equations, m: number of variables, a: augmented matrix of size n*(m+1)
const double EPS = 1e-8;

int gauss(vector<vector<double>> a, vector<double>& ans) {
  int n = (int)a.size(), m = (int)a[0].size() - 1;
  vector<int> pos(m, -1);
  double det = 1;
  int rank = 0;
  for (int r = 0, c = 0; r < n && c < m; ++c) {
    int k = r;
    for (int i = r; i < n; i++)
      if (abs(a[i][c]) > abs(a[k][c])) k = i;
    if (abs(a[k][c]) < EPS) {
      det = 0;
      continue;
    }
    swap(a[r], a[k]);
    if (r != k) det = -det;
    det *= a[r][c];
    pos[c] = r;
    for (int i = 0; i < n; i++) {
      if (i != r && abs(a[i][c]) > EPS) {
        double t = a[i][c] / a[r][c];
        for (int j = c; j <= m; j++) a[i][j] -= a[r][j] * t;
      }
    }
    ++r;
    ++rank;
  }
  ans.assign(m, 0);
  for (int i = 0; i < m; i++) {
    if (pos[i] != -1) ans[i] = a[pos[i]][m] / a[pos[i]][i];
  }
  for (int i = 0; i < n; i++) {
    double sum = 0;
    for (int j = 0; j < m; j++) sum += ans[j] * a[i][j];
    if (abs(sum - a[i][m]) > EPS) return -1;  // no solution
  }
  for (int i = 0; i < m; i++)
    if (pos[i] == -1) return 2;  // infinite solutions
  return 1;  // unique solution
}
```

+ XOR Equations

```cpp
const int N = 2010;

int gauss(int n, int m, vector<bitset<N>> a, bitset<N>& ans) {
  vector<int> pos(m, -1);
  for (int r = 0, c = 0; r < n && c < m; ++c) {
    int k = r;
    while (k < n && !a[k][c]) k++;
    if (k == n) continue;
    swap(a[r], a[k]);
    pos[c] = r;
    for (int i = 0; i < n; i++)
      if (i != r && a[i][c]) a[i] ^= a[r];
    ++r;
  }
  ans.reset();
  for (int i = 0; i < m; i++) {
    if (pos[i] != -1) ans[i] = a[pos[i]][m];
  }
  for (int i = 0; i < n; i++)
    if (((ans & a[i]).count() & 1) != a[i][m]) return -1;  // no solution
  for (int i = 0; i < m; i++)
    if (pos[i] == -1) return 2;  // infinite solutions
  return 1;  // unique solution
}
```

### Linear Basis

```cpp
struct Basis {
  vector<i64> a;
  void insert(i64 x) {
    x = minxor(x);
    if (!x) return;
    for (i64& i : a)
      if ((i ^ x) < i) i ^= x;
    a.push_back(x);
    sort(a.begin(), a.end());
  }
  bool can(i64 x) { return !minxor(x); }
  i64 maxxor(i64 x = 0) {
    for (i64 i : a) x = max(x, x ^ i);
    return x;
  }
  i64 minxor(i64 x = 0) {
    for (i64 i : a) x = min(x, x ^ i);
    return x;
  }
  i64 kth(i64 k) {  // 1st is 0
    int sz = a.size();
    if (k > (1LL << sz)) return -1;
    k--;
    i64 ans = 0;
    for (int i = 0; i < sz; i++)
      if (k >> i & 1) ans ^= a[i];
    return ans;
  }
};
```

### Chinese Remainder Theorem

```cpp
// Prerequisite: exgcd
i64 excrt(vector<i64>& m, vector<i64>& r) {
  i64 M = m[0], R = r[0];
  for (int i = 1; i < m.size(); i++) {
    auto [d, x, y] = exgcd(M, m[i]);
    if ((r[i] - R) % d) return -1;
    x = mul(x, (r[i] - R) / d, m[i] / d);
    R += x * M;
    M = M / d * m[i];
    R %= M;
  }
  return R >= 0 ? R : R + M;
}
```

### Primitive Root

```cpp
// Prerequisite: Prime Factorization
i64 primitive_root(i64 p) {
  vector<pair<i64, int>> facs = getf(p - 1);
  for (i64 i = 2; i < p; i++) {
    bool flag = true;
    for (auto& [x, _] : facs) {
      if (qk(i, (p - 1) / x, p) == 1) {
        flag = false;
        break;
      }
    }
    if (flag) return i;
  }
  return -1;
}
```

### Baby Step Giant Step

```cpp
// a^x = b (mod p), where p is a prime number
i64 BSGS(i64 a, i64 b, i64 p) {
  a %= p;
  if (!a && !b) return 1;
  if (!a) return -1;
  map<i64, i64> mp;
  i64 m = ceil(sqrt(p)), v = 1;
  for (int i = 1; i <= m; i++) {
    (v *= a) %= p;
    mp[v * b % p] = i;
  }
  i64 vv = v;
  for (int i = 1; i <= m; i++) {
    auto it = mp.find(vv);
    if (it != mp.end()) return i * m - it->second;
    (vv *= v) %= p;
  }
  return -1;
}

// Modulus can be non-prime
i64 exBSGS(i64 a, i64 b, i64 p) {
  a %= p; b %= p;
  if (a == 0) return b > 1 ? -1 : (b == 0 && p != 1);
  i64 c = 0, q = 1;
  for (;;) {
    i64 g = gcd(a, p);
    if (g == 1) break;
    if (b == 1) return c;
    if (b % g) return -1;
    ++c; b /= g; p /= g; q = a / g * q % p;
  }
  map<i64, i64> mp;
  i64 m = ceil(sqrt(p)), v = 1;
  for (int i = 1; i <= m; i++) {
    (v *= a) %= p;
    mp[v * b % p] = i;
  }
  for (int i = 1; i <= m; i++) {
    (q *= v) %= p;
    auto it = mp.find(q);
    if (it != mp.end()) return i * m - it->second + c;
  }
  return -1;
}

// Given x, b, p, find a
i64 SGSB(i64 x, i64 b, i64 p) {
  i64 g = primitive_root(p);
  return qk(g, BSGS(qk(g, x, p), b, p), p);
}
```

### Number Theory Block

```cpp
// When floor(n/i) = v, the range of i is [l,r]
vector<array<i64, 3>> quotients(i64 n) {
  vector<array<i64, 3>> res;
  i64 h = sqrt(n);
  res.reserve(2 * h - (h == n / h));
  for (i64 l = 1, v, r; l <= n; l = r + 1) {
    v = n / l;
    r = n / v;
    res.push_back({l, r, v});
  }
  return res;
}
```

### Quadratic Residue

```cpp
i64 Quadratic_residue(i64 a) {
  if (a == 0) return 0;
  i64 b;
  do b = rng() % P;
  while (qk(b, (P - 1) >> 1, P) != P - 1);
  i64 s = P - 1, t = 0, f = 1;
  while (!(s & 1)) s >>= 1, t++, f <<= 1;
  t--, f >>= 1;
  i64 x = qk(a, (s + 1) >> 1, P), inv_a = qk(a, P - 2, P);
  while (t) {
    f >>= 1;
    if (qk(inv_a * x % P * x % P, f, P) != 1) {
      (x *= qk(b, s, P)) %= P;
    }
    t--, s <<= 1;
  }
  if (x * x % P != a) return -1;
  return min(x, P - x);
}
```

### FFT & NTT & FWT

+ FFT

```cpp
const double PI = acos(-1);
using cp = complex<double>;

int n1, n2, n, k, rev[N];

void fft(vector<cp>& a, int p) {
  for (int i = 0; i < n; i++) if (i < rev[i]) swap(a[i], a[rev[i]]);
  for (int h = 1; h < n; h <<= 1) {
    cp wn(cos(PI / h), p * sin(PI / h));
    for (int i = 0; i < n; i += (h << 1)) {
      cp w(1, 0);
      for (int j = 0; j < h; j++, w *= wn) {
        cp x = a[i + j], y = w * a[i + j + h];
        a[i + j] = x + y, a[i + j + h] = x - y;
      }
    }
  }
  if (p == -1) for (int i = 0; i < n; i++) a[i] /= n;
}

void go(vector<cp>& a, vector<cp>& b) {
  n = 1, k = 0;
  while (n <= n1 + n2) n <<= 1, k++;
  a.resize(n); b.resize(n);
  for (int i = 0; i < n; i++) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (k - 1));
  fft(a, 1); fft(b, 1);
  for (int i = 0; i < n; i++) a[i] *= b[i];
  fft(a, -1);
}
```

+ NTT

```cpp
const int P = 998244353, G = 3, IG = 332748118;

int n1, n2, n, k, rev[N];

void ntt(vector<i64>& a, int p) {
  for (int i = 0; i < n; i++) if (i < rev[i]) swap(a[i], a[rev[i]]);
  for (int h = 1; h < n; h <<= 1) {
    i64 wn = qk(p == 1 ? G : IG, (P - 1) / (h << 1), P);
    for (int i = 0; i < n; i += (h << 1)) {
      i64 w = 1;
      for (int j = 0; j < h; j++, (w *= wn) %= P) {
        i64 x = a[i + j], y = w * a[i + j + h] % P;
        a[i + j] = (x + y) % P, a[i + j + h] = (x - y + P) % P;
      }
    }
  }
  if (p == -1) {
    i64 ninv = qk(n, P - 2, P);
    for (int i = 0; i < n; i++) (a[i] *= ninv) %= P;
  }
}

void go(vector<i64>& a, vector<i64>& b) {
  n = 1, k = 0;
  while (n <= n1 + n2) n <<= 1, k++;
  a.resize(n); b.resize(n);
  for (int i = 0; i < n; i++) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (k - 1));
  ntt(a, 1); ntt(b, 1);
  for (int i = 0; i < n; i++) (a[i] *= b[i]) %= P;
  ntt(a, -1);
}
```

+ FWT

```cpp
void AND(i64& a, i64& b) { a += b; }
void rAND(i64& a, i64& b) { a -= b; }

void OR(i64& a, i64& b) { b += a; }
void rOR(i64& a, i64& b) { b -= a; }

void XOR(i64& a, i64& b) {
  i64 x = a, y = b;
  a = (x + y) % P;
  b = (x - y + P) % P;
}
void rXOR(i64& a, i64& b) {
  static i64 inv2 = (P + 1) / 2;
  i64 x = a, y = b;
  a = (x + y) * inv2 % P;
  b = (x - y + P) * inv2 % P;
}

template<class T>
void fwt(vector<i64>& a, int n, T f) {
  for (int d = 1; d < n; d <<= 1) {
    for (int i = 0; i < n; i += (d << 1)) {
      for (int j = 0; j < d; j++) {
        f(a[i + j], a[i + j + d]);
      }
    }
  }
}
```

### Simpson's Numerical Integration

```cpp
double simpson(double l, double r) {
  double c = (l + r) / 2;
  return (f(l) + 4 * f(c) + f(r)) * (r - l) / 6;
}

double asr(double l, double r, double eps, double S) {
  double mid = (l + r) / 2;
  double L = simpson(l, mid), R = simpson(mid, r);
  if (fabs(L + R - S) < 15 * eps) return L + R + (L + R - S) / 15;
  return asr(l, mid, eps / 2, L) + asr(mid, r, eps / 2, R);
}

double asr(double l, double r) { return asr(l, r, EPS, simpson(l, r)); }
```

### Berlekamp–Massey

```cpp
// finding the k-th term is O(n^2 log k)
namespace BerlekampMassey {
  using V = vector<i64>;

  void up(i64 & a, i64 b) { (a += b) %= P; }

  V mul(const V& a, const V& b, const V& m, int k) {
    V r(2 * k - 1);
    for (int i = 0; i < k; i++)
      for (int j = 0; j < k; j++)
        up(r[i + j], a[i] * b[j]);
    for (int i = k - 2; i >= 0; i--) {
      for (int j = 0; j < k; j++)
        up(r[i + j], r[i + k] * m[j]);
      r.pop_back();
    }
    return r;
  }

  V pow(i64 n, const V& m) {
    int k = (int)m.size() - 1;
    assert(m[k] == -1 || m[k] == P - 1);
    V r(k), x(k);
    r[0] = x[1] = 1;
    for (; n; n >>= 1, x = mul(x, x, m, k))
      if (n & 1) r = mul(x, r, m, k);
    return r;
  }

  i64 go(const V& a, const V& x, i64 n) {
    // a: (-1, a1, a2, ..., ak).reverse
    // x: x1, x2, ..., xk
    // x[n] = sum[a[i]*x[n-i],{i,1,k}]
    int k = (int)a.size() - 1;
    if (n <= k) return x[n - 1];
    if (a.size() == 2) return x[0] * qk(a[0], n - 1, P) % P;
    V r = pow(n - 1, a);
    i64 ans = 0;
    for (int i = 0; i < k; i++) up(ans, r[i] * x[i]);
    return (ans + P) % P;
  }

  V BM(const V& x) {
    V C{-1}, B{-1};
    i64 L = 0, m = 1, b = 1;
    for (int n = 0; n < (int)x.size(); n++) {
      i64 d = 0;
      for (int i = 0; i <= L; i++) up(d, C[i] * x[n - i]);
      if (d == 0) { ++m; continue; }
      V T = C;
      i64 c = P - d * inv(b, P) % P;
      C.resize(max(C.size(), size_t(B.size() + m)));
      for (int i = 0; i < (int)B.size(); i++) up(C[i + m], c * B[i]);
      if (2 * L > n) { ++m; continue; }
      L = n + 1 - L; B.swap(T); b = d; m = 1;
    }
    reverse(C.begin(), C.end());
    return C;
  }
}
```

### Lagrange Interpolation

```cpp
// find the coefficients of f(x), O(n^2)
template <class T>
vector<T> La(vector<T> x, vector<T> y) {
  int n = x.size();
  vector<T> ret(n), sum(n);
  ret[0] = y[0], sum[0] = 1;
  for (int i = 1; i < n; i++) {
    for (int j = n - 1; j >= i; j--) {
      y[j] = (y[j] - y[j - 1]) / (x[j] - x[j - i]);
    }
    for (int j = i; j; j--) {
      sum[j] = -sum[j] * x[i - 1] + sum[j - 1];
      ret[j] += sum[j] * y[i];
    }
    sum[0] = -sum[0] * x[i - 1];
    ret[0] += sum[0] * y[i];
  }
  return ret;
}
```
