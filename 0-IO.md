## 输入 & 输出

### 特殊格式

```cpp
long double %Lf
unsigned int %u
unsigned long long %llu

cout << fixed << setprecision(15);
```

### 文件和流同步

```cpp
freopen("in.txt", "r", stdin);

ios::sync_with_stdio(false);
cin.tie(0);
```

### 程序计时

```cpp
(double)clock() / CLOCKS_PER_SEC
```

### 整行读入

```cpp
scanf("%[^\n]", s)  // 需测试是否可用
getline(cin, s)
```

### 读到文件尾

```cpp
while (cin) {}
while (~scanf) {}
```

### int128

```cpp
// 需测试是否可用
inline __int128 get128() {
  __int128 x = 0, sgn = 1;
  char c = getchar();
  for (; c < '0' || c > '9'; c = getchar()) if (c == '-') sgn = -1;
  for (; c >= '0' && c <= '9'; c = getchar()) x = x * 10 + c - '0';
  return sgn * x;
}

inline void print128(__int128 x) {
  if (x < 0) {
    putchar('-');
    x = -x;
  }
  if (x >= 10) print128(x / 10);
  putchar(x % 10 + '0');
}
```

### 读入挂

```cpp
class Scanner {
#ifdef qdd
  static constexpr int BUF_SIZE = 1;
#else
  static constexpr int BUF_SIZE = 1048576; // 1MB
#endif

  char buf[BUF_SIZE], *p1 = buf, *p2 = buf;

  char nc() {
    if (p1 == p2) {
      p1 = buf; p2 = buf + fread(buf, 1, BUF_SIZE, stdin);
      // assert(p1 != p2);
    }
    return *p1++;
  }

public:
  Scanner& operator>>(string& s) {
    s.clear();
    char c = nc();
    while (c <= 32) c = nc();
    for (; c > 32; c = nc()) s += c;
    return *this;
  }

  Scanner& operator>>(int& x) {
    x = 0;
    int sgn = 1;
    char c = nc();
    for (; c < '0' || c > '9'; c = nc()) if (c == '-') sgn = -1;
    for (; c >= '0' && c <= '9'; c = nc()) x = x * 10 + (c - '0');
    return x *= sgn, *this;
  }

  Scanner& operator>>(double& x) {
    x = 0;
    double base = 0.1;
    int sgn = 1;
    char c = nc();
    for (; c < '0' || c > '9'; c = nc()) if (c == '-') sgn = -1;
    for (; c >= '0' && c <= '9'; c = nc()) x = x * 10 + (c - '0');
    for (; c < '0' || c > '9'; c = nc()) if (c != '.') return x *= sgn, *this;
    for (; c >= '0' && c <= '9'; c = nc()) x += base * (c - '0'), base *= 0.1;
    return x *= sgn, *this;
  }
} in;
```
