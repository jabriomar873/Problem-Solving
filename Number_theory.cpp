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
ll fact[MAXN];
ll inv[MAXN];


///There are n children and m apples that will be distributed to them. 
////Your task is to count the number of ways this can be done.
//// ans=nCr(n+m-1,n-1)



ll mul(ll u ,ll v){
    return ((u%MOD)*(v%MOD))%MOD;
}
ll add(ll u, ll v){
    return ((u%MOD)+(v%MOD)+MOD)%MOD;
}
ll sub(ll u, ll v){
    return ((u%MOD)-(v%MOD)+MOD)%MOD;
}
ll power(ll a,ll b){ 
    ll res=1; 
    while(b){ 
         if(b&1) res=mul(res,a);
         a=mul(a,a); 
         b>>=1; 
    } 
    return res;
}

ll dgsum(ll x) 
{
    ll sum = 0;
    while (x > 0) {
        sum += x % 10;
        x /= 10;
    }
    return sum;
}

ll pgcd(ll a, ll b) { return b ? pgcd(b, a % b) : a;}

ll ppcm(ll a, ll b) { return (a * b) / pgcd(a, b); }



ll POW(ll a,ll b){
    ll res=1;
    while(b){
        if(b&1) res=(res*a);
        a=(a*a);
        b>>=1;
    }
    return res;
}


void precompute() {
    fact[0] = 1;
    for (ll i = 1; i < MAXN; ++i)
        fact[i] = (fact[i - 1] * 1LL * i) % MOD;
 
    inv[MAXN - 1] = power(fact[MAXN - 1], MOD - 2);
    for (ll i = MAXN - 2; i >= 0; --i)
        inv[i] = (inv[i + 1] * 1LL * (i + 1)) % MOD;
}
ll nCr(ll n, ll r) {
    if (r > n) return 0;
    ll res = (fact[n] * 1LL * inv[r]) % MOD;
    res = (res * 1LL * inv[n - r]) % MOD;
    return res;
}

ll C(ll n, ll r){
    ll p = 1, k = 1;
    if (n - r < r)r = n - r;
    if (r != 0) {
            while (r) {
                p *= n;k *= r;
                ll m = __gcd(p, k);
                p /= m;k /= m;n--;r--;
                }
            }
    else p = 1;
    return p;
}
 
vll DD[MAXN];
ll SPF[NN];
void sieve() {
    for (ll i = 0; i <NN; ++i) SPF[i] = i;
    for (ll i = 2; i <NN; ++i) {
        if (SPF[i] != i)continue;
        for (ll j = 2 * i; j <NN; j += i)
            SPF[j] = min(SPF[j], i);
    }
}
mpll factorize(ll x) {
    mpll res;
    while (x > 1) {
        res[SPF[x]]++;
        x /= SPF[x];
    } return res;
}

ll nbr_divieur(ll n)
{
    mpll mp;
    while (n%2== 0) mp[2]++, n = n/2;
    for (ll i = 3; i*i <= n; i = i+2)
    {
        while (n%i == 0)
        {   
            mp[i]++;
            n = n/i;
        }
    }
    if (n > 2) mp[n]++;
    ll nbr=1;
    for(auto m:mp) nbr*=(m.second+1);
    return nbr;
}
ll DA[NN];
void nbr_Divisors(){
    for(ll i=1;i<=NN;i++)
        for(ll j=i;j<=NN;j+=i)   ///also we can use to find diviseurs en n*log(log(n)) juste divs[j].pb(i);
            DA[j]++;
}
ll gcd_extended(ll a,ll b,ll &x,ll &y)  /// return PGCD and the solution (x,y) for (E):a*x+b*y+c;
{
    if(a==0) 
    {
        x=0;                     /// !!! we have a solution if c%gcd==0
        y=1;
        return b;
    }
    ll x1,y1;
    ll gcd=gcd_extended(b%a,a,x1,y1);
    x=y1-(b/a)*x1;            /// solve a*x+b*y+c=0 : Bezout equation  !!!! +c
    y=x1;                      /// x*=-c/gcd and y*=-c/gcd  because we are solving a*x+b*y=gcd(a,b)
    return gcd;                /// when we multiple by c/gcd(a,b) we get a*x+b*y=c
}

long long  phi(long long n)
{
	long long res = n;
	for (long long p = 2; p * p <= n; ++p)
	{
		if (n % p == 0)
			res -= res / p;
		while (n % p == 0)
			n /= p;
	}
	if (n > 1)
		res -= res / n;
	return res;
}

ll inclusion_exclusion(vll&v,ll n,ll k){
    ll odd=0;
    ll even=0;
    for(ll i=1;i<(1<<k);i++){
        ll p=1;
        bool test=true;
        for(ll j=0;j<k;j++){
            if((1<<j)&i) {
                ll r=INF/p;
                if(v[j]>r){
                    test=false;
                    break;
                }
                p*=v[j];
            }
        }
        if(!test) continue;
        ll nbr=__builtin_popcountll(i);
        if(nbr%2) odd+=(n/p);
        else even+=(n/p);
    }
    return odd-even;
}

const ll N=2e6+2;
ll mob[N];
bool prime[N];
vll D[N];
void moebius(){
    fill(mob, mob + N, 1);
    memset(prime + 2, 1, sizeof(prime) - 2);
    mob[0] = 0;
    mob[2] = -1;
    for (ll i = 4; i < N; i += 2){
        mob[i] *= (i & 3) ? -1 : 0;
        prime[i] = 0;
    }
    for (ll i = 3; i < N; i += 2){
        if (prime[i]){
            mob[i] = -1;
            for (ll j = 2 * i; j < N; j += i){
                mob[j] *= j % (1LL * i * i) ? -1 : 0;
                prime[j] = 0;
            }
        }
    }
    for(ll i = 2; i < N; i++){
        if(mob[i] == 0) continue;
        for(ll j = i; j < N; j += i){
            D[j].push_back(i);
        }
    }
}
ll freq[N];
//all solution of ax+by=c in [minx,maxx], [miny,maxy]
int gcd(int a, int b, int &x, int &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int x1, y1;
    int d = gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}

bool solve_any(int a, int b, int c, int &x, int &y, int &g) {
    if (a == 0 && b == 0) {
        g = x = y = 0;
        return true;
    }
    g = gcd(abs(a), abs(b), x, y);
    if (c % g != 0) {
        return false;
    }
    x *= (c / g);
    y *= (c / g);
    if (a < 0) x = -x;
    if (b < 0) y = -y;
    return true;
}

void shift(int &x, int &y, int a, int b, int cnt) {
    x += cnt * b;
    y -= cnt * a;
}

int solve_all(int a, int b, int c, int minx, int maxx, int miny, int maxy) {
    int x, y, g;
    if (!solve_any(a, b, c, x, y, g))
        return 0;
    a /= g;
    b /= g;
    int sign_a = a > 0 ? +1 : -1;
    int sign_b = b > 0 ? +1 : -1;
    shift(x, y, a, b, (minx - x) / b);
    if (x < minx) {
        shift(x, y, a, b, sign_b);
    }
    if (x > maxx)
        return 0;
    int lx1 = x;
    shift(x, y, a, b, (maxx - x) / b);
    if (x > maxx) {
        shift(x, y, a, b, -sign_b);
    }
    int rx1 = x;
    shift(x, y, a, b, (y - miny) / a);
    if (y < miny) {
        shift(x, y, a, b, -sign_a);
    }
    if (y > maxy)
        return 0;
    int lx2 = x;

    shift(x, y, a, b, (y - maxy) / a);
    if (y > maxy) {
        shift(x, y, a, b, sign_a);
    }
    int rx2 = x;
    if (lx2 > rx2)
        swap(lx2, rx2);
    int lx = max(lx1, lx2);
    int rx = min(rx1, rx2);
    if (lx > rx) return 0;
    return (rx - lx) / abs(b) + 1;
}

void solve() {
    int a, b, c;
    cin >> a >> b >> c;
    int minx, miny, maxx, maxy;
    cin >> minx >> maxx >> miny >> maxy;
    if (maxx < minx || maxy < miny) {
        cout << 0 << endl;
        exit(0);
    }
    if (a == 0 && b == 0) {
        if (c == 0) {
            cout << (maxy - miny + 1) * (maxx - minx + 1) << endl;
        } else {
            cout << 0 << endl;
        }
    } else if (a == 0 || b == 0) {
        if (a == 0) {
            if (c % b == 0 && (c / b) >= miny && (c / b) <= maxy) {
                cout << maxx - minx + 1 << endl;
            } else {
                cout << 0 << endl;
            }
        } else {
            if (c % a == 0 && (c / a) >= minx && (c / a) <= maxx) {
                cout << maxy - miny + 1 << endl;
            } else {
                cout << 0 << endl;
            }
        }
    } else
        cout << solve_all(a, b, c, minx, maxx, miny, maxy) << endl;
}
