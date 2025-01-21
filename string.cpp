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

bool is_palindrome(const string &s, int start, int len) {
    int end = start + len - 1;
    while (start < end) {
        if (s[start] != s[end]) return false;
        start++;
        end--;
    }
    return true;
}

/*-------hashing----------*/
struct custom_hash {
    static uint64_t splitmix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }
    template <class T>
    size_t operator()(const T& x) const {
        return splitmix64(hash<T>()(x) + 0x9e3779b9);
    }
    template <class T, class U>
    size_t operator()(const pair<T, U>& p) const {
        auto hash1 = hash<T>()(p.first);
        auto hash2 = hash<U>()(p.second);
        return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
    }
    template <class T>
    struct hash {
        size_t operator()(const T& x) const {
            return std::hash<T>{}(x);
        }
    };
};
struct Hashing { //// reversible but you can not get the hash of a substring s[i,j]
private:
    int mod1 = 1e9 + 7, mod2 = 2e9 + 11;
    ll base1, base2, h1, h2, inv1, inv2, *pw1, *pw2, len;
    deque<char> d;
    ll power(ll a, ll b, ll m) {
        ll ans = 1;
        while (b > 0) {
            if (b & 1) ans = (ans * a) % m;
            a = (a * a) % m;
            b >>= 1;
        }
        return ans;
    }
public:
    Hashing(ll sz, ll x = 31, ll y = 37) {
        base1 = x;
        base2 = y;
        h1 = h2 = len = 0;
        inv1 = power(x, mod1 - 2, mod1);
        inv2 = power(y, mod2 - 2, mod2);
        pw1 = new ll[sz + 1];
        pw2 = new ll[sz + 1];
        pw1[0] = pw2[0] = 1;
        for (int i = 1; i <= sz; i++) {
            pw1[i] = (x * pw1[i - 1]) % mod1;
            pw2[i] = (y * pw2[i - 1]) % mod2;
        }
    }
    void push_back(char x) {
        x = x - 'a' + 1;
        h1 = (h1 * base1) % mod1;
        h1 = (h1 + x) % mod1;
        h2 = (h2 * base2) % mod2;
        h2 = (h2 + x) % mod2;
        len++;
        d.emplace_back(x);
    }
    void push_front(char x) {
        x = x - 'a' + 1;
        h1 = (h1 + (x * pw1[len]) % mod1) % mod1;
        h2 = (h2 + (x * pw2[len]) % mod2) % mod2;
        len++;
        d.emplace_front(x);
    }
    void pop_back() {
        if (len == 0) return;
        char x = d.back();
        d.pop_back();
        h1 = (h1 - x + mod1) % mod1;
        h1 = (h1 * inv1) % mod1;
        h2 = (h2 - x + mod2) % mod2;
        h2 = (h2 * inv2) % mod2;
        len--;
    }
    void pop_front() {
        if (len == 0) return;
        char x = d.front();
        d.pop_front();
        len--;
        h1 = ((h1 - x * pw1[len] % mod1) + mod1) % mod1;
        h2 = ((h2 - x * pw2[len] % mod2) + mod2) % mod2;
    }
    void clear() {
        h1 = h2 = len = 0;
        d.clear();
    }
    bool operator==(const Hashing &H) const {
        return H.h1 == h1 && H.h2 == h2;
    }
    string GetString() {
        return string(d.begin(), d.end());
    }
    pair<ll,ll> GetHash() {
        return {h1, h2};
    }
};


struct Hashing { /// not reversible but you can  get the hash of a substring s[i,j]
private:
    ll mod1 = 1e9 + 7, mod2 = 2e9 + 11;
    ll base1, base2, inv1, inv2, len;
    vector<ll> h1, h2, pw1, pw2;
    ll power(ll a, ll b, ll m) {
        ll ans = 1;
        while (b > 0) {
            if (b & 1) ans = (ans * a) % m;
            a = (a * a) % m;
            b >>= 1;
        }
        return ans;
    }
public:
    Hashing(ll sz, ll x = 31, ll y = 37) {
        base1 = x;
        base2 = y;
        len = 0;
        inv1 = power(x, mod1 - 2, mod1);
        inv2 = power(y, mod2 - 2, mod2);
        h1.resize(sz + 1);
        h2.resize(sz + 1);
        pw1.resize(sz + 1);
        pw2.resize(sz + 1);
        pw1[0] = pw2[0] = 1;
        for (int i = 1; i <= sz; i++) {
            pw1[i] = (pw1[i - 1] * x) % mod1;
            pw2[i] = (pw2[i - 1] * y) % mod2;
        }
    }
    void addChar(char x) {
        x = x - 'a' + 1;
        h1[len + 1] = (h1[len] * base1 + x) % mod1;
        h2[len + 1] = (h2[len] * base2 + x) % mod2;
        len++;
    }
    pll substring_hash(ll i, ll j) {
        ll hash_value_1 = (h1[j + 1] - (h1[i] * pw1[j - i + 1]) % mod1 + mod1) % mod1;
        ll hash_value_2 = (h2[j + 1] - (h2[i] * pw2[j - i + 1]) % mod2 + mod2) % mod2;
        return {hash_value_1, hash_value_2};
    }
    bool is_palindrome(ll i, ll j) {
        return substring_hash(i, j) == substring_hash(len - j - 1, len - i - 1);
    }
    void clear() {
        len = 0;
        fill(h1.begin(), h1.end(), 0);
        fill(h2.begin(), h2.end(), 0);
    }
};

/// Haashing H(s.length);
//// usse an ordred_set for each pair of hash 
//// synx : ordred_set<pll,pair_hash> st; access in O(1)


/*-----------substring----------*/
bool isSubstring(const string& str, const string& substr) 
{
    size_t found = str.find(substr);
    return found != string::npos;
}
bool compareStrings(const string& a, const string& b) 
{
    return a.size() < b.size() && isSubstring(b, a);
}

/**KMP**/
void computeLPSArray(string pat,ll M,ll lps[]){
    ll len = 0,i=1;
    lps[0] = 0;   /// lps[i]=k : s[0..k-1]=s[i-k+1..i]
    while (i < M){
            if (pat[i] == pat[len]){len++;lps[i] = len;i++;}
            else {if (len != 0){len = lps[len - 1];}
            else {lps[i] = len;i++;}}}}
 
ll KMP(string pat, string txt){
    ll M = pat.length(),N = txt.length(),lps[M],j=0,i=0,res=0;
    computeLPSArray(pat, M, lps);
    while (i < N){
            if (pat[j] == txt[i]){j++;i++;}
            if (j == M){j = lps[j - 1];res++;}
            else if (i < N && pat[j] != txt[i]){
                        if (j != 0)j = lps[j - 1];
                        else i = i + 1;
                    }
    }
    return res;  //// return freq of pat in txt
}



ll countPalindromicSubstrings(const string& s) {  /// Manacher'alogorithme
    // Step 1: Transform the string s into T
    string T = "#";
    for (char c : s) {
        T += c;
        T += "#";
    }
    ll n = T.length();
    // Step 2: Initialize variables
    vll P(n, 0);  // Array to store the radius of palindromes
    ll C = 0, R = 0;     // Current center and right edge
    ll totalPalindromes = 0;
    // Step 3: Expand around each center
    for (ll i = 0; i < n; ++i) {
        ll mirror = 2 * C - i; // Mirror of i
        // Check if i is within the current palindrome
        if (i < R)
            P[i] = min(R - i, P[mirror]);
        // Expand around the center i
        while (i + P[i] + 1 < n && i - P[i] - 1 >= 0 && T[i + P[i] + 1] == T[i - P[i] - 1])
            P[i]++;
        // Update center and right edge
        if (i + P[i] > R) {
            C = i;
            R = i + P[i];
        }
        // Count palindromes in the original string
        totalPalindromes += (P[i] + 1) / 2;
    }
    return totalPalindromes;
}

struct Manacher {
    // Function to pad the string with a delimiter (default is '#')
    string getPaddedString(const string& s, string delim = "#") {
        string paddedString = "#";
        for (const char& c : s) {
            paddedString += c;
            paddedString += delim;
        }
        return paddedString;
    }
    // Function to find the position of the longest palindrome
    long getMaxPos(const vll& span) {
        ll maxlen = 1, pos = 1, n = span.size();
        for (ll i = 2; i < n; i++) {
            if (span[i] > maxlen) {
                maxlen = span[i];
                pos = i;
            }
        }
        return pos;
    }
    // Function to get the unpadded string (longest palindrome)
    string getUnpaddedString(const string& s, const ll& maxPos, const ll& maxlen) {
        string ans;
        ll start = maxPos - maxlen + 1, end = maxPos + maxlen - 1;
        for (ll i = start; i <= end; i++) {
            if (s[i] != '#') ans += s[i];
        }
        return ans;
    }
    // Main Manacher's algorithm function
    string manacher(const string& s) {
        string str = getPaddedString(s);
        ll n = str.length();
        vll span(n, 1);
        ll c = 0, r = 0;
        // Manacher's algorithm
        for (ll i = 1; i < n; i++) {
            if (i <= r && 2 * c - i >= 0) {
                span[i] = min(r - i + 1, span[2 * c - i]);
            }
            while (i - span[i] >= 0 && i + span[i] < n && str[i - span[i]] == str[i + span[i]]) {
                span[i]++;
            }
            if (i + span[i] - 1 > r) {
                c = i;
                r = i + span[i] - 1;
            }
        }
        ll maxPos = getMaxPos(span);
        ll maxlen = span[maxPos];
        return getUnpaddedString(str, maxPos, maxlen);
    }
};
vector<ll> z_function(const string& s) {
	ll n = s.length();
	vector<ll> z(n);

	for (ll i = 1, l = 0, r = 0; i < n; i++) {
		if (i <= r) {
			z[i] = min(z[i-l], r-i+1);
		}
		while (i + z[i] < n and s[z[i]] == s[i+z[i]]) { ///role=Finding the longest substring which is also a prefix:
			++z[i];
		}
		if (i + z[i] - 1 > r) {
			l = i, r = i + z[i] - 1;
		}
	}
	return z;    // :  z[i]=k it means s[i,....,i+k-1]=s[0,...,k-1] ,z[0]=n trivial value
                 // we can use z_function like KMP or more better 
                 // i can return all the pos of patern using z_function
}

//applications are same as prefix function