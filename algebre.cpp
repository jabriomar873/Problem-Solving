 #include <bits/stdc++.h>
 using namespace std;
 #define Try ios::sync_with_stdio(0); cin.tie(0); cout.tie(0); cout << std::setiosflags(std::ios::fixed) << std::setprecision(20);
 #define dd long double
 #define endl "\n"
 #define ll long long
 #define ull unsigned long long
 #include <ext/pb_ds/assoc_container.hpp>
 #include <ext/pb_ds/tree_policy.hpp>
 
 using namespace __gnu_pbds;
 template <class T> using ordered_set=tree<T, null_type,less<T>, rb_tree_tag, tree_order_statistics_node_update>;
 #define all(x) x.begin(), x.end()
 #define rall(a) (a).rbegin(), (a).rend()
 #define YES cout << "YES" <<endl;
 #define NO cout << "NO" <<endl;
 #define Yes cout << "Yes" <<endl;
 #define No cout << "No" <<endl;
 #define lire(v) for (auto &u : v) cin >> u;
 #define aff(v) for (auto u : v) cout << u << " "; cout <<endl;
 #define Unique(v) v.erase(unique(all(v)), v.end())
 #define sq(a) (a)*(a)
 #define ff first
 #define ss second
 #define vll vector<ll>
 #define vpll vector<pair<ll, ll>>
 #define vbo vector<bool>
 #define mpll map<ll,ll>
 #define sll set<ll>
 #define pll pair<ll, ll>
 #define vvll vector<vll>
 
 
const double PI = acos(-1);
const dd EPS = 1e-9;
const ll MOD = 1e9 +7;
///const ll MOD = 998244353;
const long long INF =1e18;
const ll MAXN =2e5+2;
const ll NN=1e6+1;



vector<vector<char>> rotation(vector<vector<char>>&v ,ll n)
{
    vector<vector<char>> w(n,vector<char>(n));
    for(ll i=0;i<n;i++)
    {
        for(ll j=0;j<n;j++)     ///rotation de 90 dans le sens des aiguilles
        {
            w[j][n-i-1]=v[i][j];
        }
    }
    return w;
}


#define vvi vector<vector<char>>

void verticalShift(vvi &matrix,ll H,ll W,ll k) {
    k = k % H; // To avoid unnecessary full cycles
    if (k == 0) return;

    vvi temp(H, vector<char>(W));
    for (ll j = 0; j < W; j++) {
        for (ll i = 0; i < H; i++) {    //// we shifft each column from left for the vector
            temp[i][j] = matrix[(i + k) % H][j];
        }
    }
    matrix = temp;
}

void horizontalShift(vvi &matrix,ll H,ll W,ll r) {
    r = r % W; // To avoid unnecessary full cycles
    if (r == 0) return;

    vvi temp(H, vector<char>(W));
    for (ll i = 0; i < H; i++) {     //// we shifft each ligne from left for the vector
        for (ll j = 0; j < W; j++) {
            temp[i][j] = matrix[i][(j + r) % W];
        }
    }
    matrix = temp;
}

void performShifts(vvi &matrix,ll k,ll r) {
    ll H = matrix.size();
    ll W = matrix[0].size();
    verticalShift(matrix, H, W, k);
    horizontalShift(matrix, H, W, r);
}


// Function to multiply two matrices (mat1 * mat2) under modulo MOD
vvll produit(const vvll &mat1, const vvll &mat2, ll mod = MOD) {
    ll n = mat1.size();
    vvll result(n, vll(n, 0));
    for (ll i = 0; i < n; ++i) {
        for (ll j = 0; j < n; ++j) {
            for (ll k = 0; k < n; ++k) {
                result[i][j] = (result[i][j] + (mat1[i][k] % MOD * mat2[k][j] % MOD) % MOD) % MOD;
            }
        }
    }
    return result;
}
// Function to perform matrix exponentiation (mat^exp) under modulo MOD
vvll expon(vvll mat, ll exp, ll mod = MOD) {
    ll n = mat.size();
    vvll result(n, vll(n, 0));
    // Initialize result as the identity matrix
    for (ll i = 0; i < n; ++i) {
        result[i][i] = 1;
    }
    // Exponentiate the matrix using binary exponentiation
    while (exp > 0) {
        if (exp % 2 == 1) {
            result =produit(result, mat, mod);
        }
        mat =produit(mat, mat, mod);
        exp /= 2;
    }
    return result;
}

/*----FFT------*/
//recursive

using cd = complex<dd>;
const dd PI = acos(-1);

void fft(vector<cd> &a, bool inverse) {
    int n = a.size();
    if (n == 1) return;
    vector<cd> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; i++) {
        even[i] = a[2 * i];
        odd[i] = a[2 * i + 1];
    }
    fft(even, inverse);
    fft(odd, inverse);

    dd theta = 2 * PI / n * (inverse ? -1 : 1);
    cd w_j(1);
    cd w_theta(cos(theta), sin(theta));
    for (int j = 0; j < n / 2; j++) {
        a[j] = even[j] + w_j * odd[j];
        a[j + n / 2] = even[j] - w_j * odd[j];
        if (inverse) {
            a[j] /= 2;
            a[j + n / 2] /= 2;
        }
        w_j *= w_theta;
    }
}

vector<int> multiply(vector<int> &a, vector<int> &b) {
    vector<cd> fa(all(a)), fb(all(b));
    int n = 1;
    while (n < a.size() + b.size())
        n <<= 1;
    fa.resize(n);
    fb.resize(n);
    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++) {
        fa[i] *= fb[i];
    }
    fft(fa, true);
    vector<int> ans(n);
    for (int i = 0; i < n; i++)
        ans[i] = round(fa[i].real());
    return ans;
}

//iterative

using cd = complex<dd>;

void fft(vector<cd> &a, bool inverse) {
    int n = a.size();

    for (int i = 1, j = 0; i < n; i++) {
        int bit = (n >> 1);
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        dd theta = 2 * PI / len * (inverse ? -1 : 1);
        cd w_theta(cos(theta), sin(theta));
        for (int i = 0; i < n; i += len) {
            cd w_j(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i + j];
                cd v = a[i + j + len / 2] * w_j;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w_j *= w_theta;
            }
        }
    }
    if (inverse) {
        for (auto &c:a)
            c /= n;
    }
}

vector<int> multiply(vector<int> &a, vector<int> &b) {
    vector<cd> fa(all(a)), fb(all(b));
    int n = 1;
    while (n < a.size() + b.size())
        n <<= 1;
    fa.resize(n);
    fb.resize(n);
    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++) {
        fa[i] *= fb[i];
    }
    fft(fa, true);
    vector<int> ans(n);
    for (int i = 0; i < n; i++)
        ans[i] = round(fa[i].real());
    return ans;
}

//use case
//check for every len(b) in string A times ch1 and ch2 are equal (string character matching)
vector<int> find(string &a,string &b,char ch1,char ch2) {
    int n = a.size();
    int m = b.size();
    vector<int> poly1(n, 0), poly2(m, 0);
    for (int i = 0; i < n; i++) {
        poly1[i] = (a[i] == ch1);
    }
    for (int i = 0; i < m; i++) {
        poly2[m - i - 1] = (b[i] == ch2);
    }
    vector<int> ans = multiply(poly1, poly2); //fft
    vector<int> to_return;
    for (int i = m - 1; i <= n - 1; i++)
        to_return.push_back(ans[i]);
    return to_return;
}


class mat {
public:
    vector<vector<int>> m;
    int n;

    explicit mat(int _n) {
        n = _n;
        m.resize(n, vector<int>(n, 0));
    }

    mat operator*(mat const &b) {
        mat ans(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    (ans.m[i][j] += m[i][k] * b.m[k][j]) %=MOD;
                }
            }
        }
        return ans;
    }

    vector<int> operator*(vector<int>&b) const {
        vector<int> ans(n, 0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                (ans[i] += m[i][j] * b[j])%=MOD;
            }
        }
        return ans;
    }

    mat power(int p) {
        mat ans(n);
        mat curr = *this;
        for (int i = 0; i < n; i++) {
            ans.m[i][i] = 1;
        }
        while (p) {
            if (p & 1) {
                ans = ans * curr;
            }
            p /= 2;
            curr = curr * curr;
        }
        return ans;
    }
};

const int mod = 998244353;

struct Mat {
  int n, m;
  vector<vector<int>> a;
  Mat() { }
  Mat(int _n, int _m) {n = _n; m = _m; a.assign(n, vector<int>(m, 0)); }
  Mat(vector< vector<int> > v) { n = v.size(); m = n ? v[0].size() : 0; a = v; }
  inline void make_unit() {
    assert(n == m);
    for (int i = 0; i < n; i++)  {
      for (int j = 0; j < n; j++) a[i][j] = i == j;
    }
  }
  inline Mat operator + (const Mat &b) {
    assert(n == b.n && m == b.m);
    Mat ans = Mat(n, m);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < m; j++) {
        ans.a[i][j] = (a[i][j] + b.a[i][j]) % mod;
      }
    }
    return ans;
  } 
  inline Mat operator - (const Mat &b) {
    assert(n == b.n && m == b.m);
    Mat ans = Mat(n, m);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < m; j++) {
        ans.a[i][j] = (a[i][j] - b.a[i][j] + mod) % mod;
      }
    }
    return ans;
  }
  inline Mat operator * (const Mat &b) {
    assert(m == b.n);
    Mat ans = Mat(n, b.m);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < b.m; j++) {
        for(int k = 0; k < m; k++) {
          ans.a[i][j] = (ans.a[i][j] + 1LL * a[i][k] * b.a[k][j] % mod) % mod;
        }
      }
    }
    return ans;
  }
  inline Mat pow(long long k) {
    assert(n == m);
    Mat ans(n, n), t = a; ans.make_unit();
    while (k) {
      if (k & 1) ans = ans * t;
      t = t * t;
      k >>= 1;
    }
    return ans;
  }
  inline Mat& operator += (const Mat& b) { return *this = (*this) + b; }
  inline Mat& operator -= (const Mat& b) { return *this = (*this) - b; }
  inline Mat& operator *= (const Mat& b) { return *this = (*this) * b; }
  inline bool operator == (const Mat& b) { return a == b.a; }
  inline bool operator != (const Mat& b) { return a != b.a; }
};

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  int n; long long k; cin >> n >> k;
  Mat a(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      cin >> a.a[i][j];
    }
  } 
  Mat ans = a.pow(k);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      cout << ans.a[i][j] << ' ';
    }
    cout << '\n';
  }
  return 0;
}
// https://judge.yosupo.jp/problem/pow_of_matrix
