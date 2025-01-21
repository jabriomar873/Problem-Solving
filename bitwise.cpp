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



string ll__bin(ll n, ll nbr_bit) {
    string s = bitset<64>(n).to_string();
    return s.substr(64 - nbr_bit);
}

int countBits(int mask) {	//	O(bits Count)		__builtin_popcount
	int ret = 0;
	while (mask){
		mask &= (mask-1);
		++ret;	// Simply remove the last bit and so on
	}
	return ret;
}

// Get bit
// rq ind de droite vers gauche (0--->)
int getBit(int num, int idx) {
  return ((num >> idx) & 1) == 1;  // get (idx+1) ime bit from right 
}

int setBit1(int num, int idx) 
{
	return num | (1<<idx);      // set(idx+1) ime bit from right 1
}

int setBit0(int num, int idx) {          // set(idx+1) ime bit from right 0
	return num & ~(1<<idx);				// 110100, idx = 3  -->  110100 & ~000100 = 110100 & 111011
}

int flipBit(int num, int idx) {
	return num ^(1<<idx);         /// inverse le (idx+1) bit from right 
}
string conversion_bin_gray(string s) /// s string binaire naturel
{
       string ch="";
       ll n=s.size();
       reverse(s.begin(),s.end());
       for(ll i=0;i<n;i++)
       {
          if(i!=n-1)
          {
             if(s[i]==s[i+1]) ch+='0';
             else ch+='1';
             
          }
          else ch+=s[n-1];
       }
       reverse(ch.begin(),ch.end());
       return ch;
}

///  X & ~(X-1) 	= 01101001000 & 10010111000 = 00000001000	value of 1<<SmaintestBitIdx

// (A ^ B ^ C ^ D ^ E) ^ (A ^ B ^ C) = D ^ E

// X >> 1	 = 	1001000000 			Equals X/2
// X >> 2	 = 	100100000 			Equals X/4
// X << 1	 =  100100000000		Equals X*2
// X << 2	 =  1001000000000		Equals X*4
void STL() {
	const int N = 20;		// const
	string s = "000111";
	bitset<N> x(s);			// 00000000000000000111
	x.set();				// 11111111111111111111
	x.flip();				// 00000000000000000000

	x = 10;					// 00000000000000001010
	x |= 3;					// 00000000000000001011
	x = (x<<3);				// 00000000000001011000
	x = ~x;					// 11111111111110100111
	x.set(15, 0);			// 11110111111110100111
	x.set(15);				// 11111111111110100111
	x.flip(0);				// 11111111111110100110
	x.count();				// Returns the number of bits that are set.
	x.any();				// Returns true if any bits are set.
	x.none();				// Returns true if no bits are set.
	x.test(15);
	x.to_ulong();			// Returns an unsigned long represent mask

	// The most interesting
	if(x[2] == 0)
		;

	x[0] = 1;	// Set bit from most right to 1
	x[N-1] = 0;	// Set bit from most left to 0

	cout<<x<<"\n";			// display a string of N bits
}

////////////////////////////////////////
string grayCode(int i, int k) {
    // Step 1: Compute the Gray code for i
    int gray = i ^ (i >> 1);
 
    // Step 2: Convert the Gray code to a binary string of length k
    std::string binary = std::bitset<64>(gray).to_string(); // assuming k <= 64
    binary = binary.substr(64 - k); // Get the last k bits
 
    // Step 3: Replace '0' with 'W' and '1' with 'B'
    for (char &c : binary) {
        c = (c == '0') ? 'W' : 'B';
    }
 
    return binary;
}