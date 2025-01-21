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

//  #pragma GCC optimize("unroll-loops")
//  #pragma GCC optimize("unroll-all-loops")
//  #pragma GCC optimize("O3")
//  #pragma GCC optimize("Ofast")





mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
ll rnd(ll B) { return (ull)rng() % B; }

/////  -------LCS--------


vector<char> longestCommonSubsequence(const vector<char>& a, const vector<char>& b) {
    ll m = a.size();
    ll n = b.size();
    vector<vll> dp(m + 1,vll(n + 1, 0));
    // Fill the dp array
    for (ll i = 1; i <= m; ++i) {
        for (ll j = 1; j <= n; ++j) {
            if (a[i - 1] == b[j - 1]) {
                dp[i][j] = dp[i - 1][j - 1] + 1;
            } else {
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
            }
        }
    }

    // Traceback to find the LCS
    vector<char> lcs;
    ll i = m, j = n;
    while (i > 0 && j > 0) {
        if (a[i - 1] == b[j - 1]) {
            lcs.push_back(a[i - 1]);
            --i;
            --j;
        } else if (dp[i - 1][j] > dp[i][j - 1]) {
            --i;
        } else {
            --j;
        }
    }
    reverse(lcs.begin(), lcs.end());
    return lcs;
}


///// remove non necessary maks:
    /*
    vbo exist((1LL<<hh),false);
    for(auto mask:masks){
       exist[mask]=true;
    }
    for (ll i = 0; i <xx;++i) {
        for (ll subMask = (masks[i] - 1) & masks[i]; subMask > 0; subMask = (subMask - 1) & masks[i]) {
            exist[subMask] = false;
        }
    }
    */