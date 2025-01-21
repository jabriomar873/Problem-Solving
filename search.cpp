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



ll bs(vector<ll>& arr, ll target) 
{
    ll low = 0;
    ll high = arr.size() ;
    
    while (low < high) {
        ll mid = low + ((high - low) / 2);
        
        if (arr[mid]<=target)low=mid+1;  //// nbr<=     the best
        else high=mid;
    }
    return low;
}

bool bbss(vector<ll>& v,ll x,ll n)
{
    ll ft=0,sd=n-1,md;
    while (ft<=sd)
    {
        md=ft+(sd-ft)/2;
         if(v[md]==x) return true; ///// l'element existe  dans le tableau
         else if(v[md]<x) ft=md+1;
         else sd=md-1;
    }
    return false; ///// l'element n'existe pas dans le tableau
}


ll bs____intv(vector<ll>& arr, ll target,ll deb ,ll fin) 
{
    ll low =deb; /// 0
    ll high =fin+1; /// arr.size() 
    
    while (low < high) {
        ll mid = low + ((high - low) / 2);   /// also the most power!!!!!!
        
        if (arr[mid]<=target)low=mid+1;
        else high=mid;
    }
    
    return low;
}

vll subset;
ll n;
void search(ll k) {  /// generate all subset for the range [0..n-1]
    if (k == n) {
        ///aff(subset); process subset
    } else {
        search(k+1);
        subset.push_back(k);
        search(k+1);
        subset.pop_back();
    }
 } //// search(0);

/// using looping
  /*for (ll b = 0; b < (1<<n); b++) {
    vll subset;
    for (ll i = 0; i < n; i++) {
      if (b&(1<<i)) subset.push_back(i);
    }
    ///aff(subset); process subset
 }
 */
 
    /*
    vll permutation(n);
    iota(all(permutation),0); // Fill permutation with 0, 1, 2, ..., n-1
    do {
        //processus:
    } while (next_permutation(permutation.begin(), permutation.end()));
    */


// Function to recursively generate all product strings or number
void generateProducts(const vector<vector<char>>& groups, int index,string currentProduct, set<string>& products) {
        // Base case: if we have processed all groups
        if (index == groups.size()) {
            products.insert(currentProduct);
            return;
        }
        // Recursive case: iterate through all elements in the current group
        for (char element : groups[index]) {
            generateProducts(groups, index + 1, currentProduct + element, products);
        }
}

// Function to recursively generate all products
  void generateProducts(const vector<vector<int>>& groups, int index, long long currentProduct, set<long long>& products) {
    // Base case: if we have processed all groups
    if (index == groups.size()) {
        products.insert(currentProduct);
        return;
    }

    // Recursive case: iterate through all elements in the current group
    for (int element : groups[index]) {
        generateProducts(groups, index + 1, currentProduct * element, products);
    }
}

// Function to generate and print combinations of edges
void generateProducts(const vector<vpll> &groups, ll index, set<pll> &st, map<pll, ll> &remp,ll &k,ll&m) {
    if (!k) return; // If k combinations are printed, stop
    if (index == groups.size()) { // If all groups are processed
        string ch(m, '0');
        for (const pll &p : st) {
            ll ind = remp[p];
            ch[ind - 1] = '1'; // Mark the edges present in the combination
        }
        cout << ch << endl;
        k--;
        return;
    }
    for (const pll &element : groups[index]) {
        st.insert(element);
        generateProducts(groups, index + 1, st, remp,k,m); // Recurse for next group
        st.erase(element);
    }
}
    /*
    int l=0;
    int r=1e9+1;
    while(r-l>3){
        int mid=(l+r)/2;
        if(f(mid)>f(mid+1))
            l=mid+1;
        else
            r=mid;
        dbg(l,r);
    }
    int ans=inf;
    for(int i =l;i<=r;i++)
        ans=min(ans,f(i));
    cout<<ans<<endl;
    */
