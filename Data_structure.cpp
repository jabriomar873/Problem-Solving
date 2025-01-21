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


/*-----------SGT___without_lazy------*/
class segte{
public:
    struct node{
        int val;
        node(){
            val=0;
        }
    };
    int n;
    vector<node> tree;
    vector<int> a;
    node neutral;

    void init(int N) {
        n = N;
        tree.resize(4 * n + 1);
        //default values
        a.assign(n, 0);
        neutral.val=0;
    }

    void put(vector<int> &val) {
        a = val;
    }

    //merge function
    void merge(node &curr, node &left, node &right) {
        curr.val = left.val + right.val;
    }

    //for leaf
    void single(node &curr, int idx) {
        curr.val = a[idx];
    }

    void build(int index, int ss, int se) {
        if (ss == se) {
            single(tree[index], ss);
            return;
        }
        int mid = (ss + se) / 2;
        build(2 * index, ss, mid);
        build(2 * index + 1, mid + 1, se);
        merge(tree[index], tree[2 * index], tree[2 * index + 1]);
    }

    void build() {
        build(1, 0, n - 1);
    }

    node query(int index, int ss, int se, int qs, int qe) {
        if (qs > se || qe < ss) return neutral;
        if (qs <= ss && qe >= se) return tree[index];
        int mid = (ss + se) / 2;
        node left = query(2 * index, ss, mid, qs, qe);
        node right = query(2 * index + 1, mid + 1, se, qs, qe);
        node mer;
        merge(mer, left, right);
        return mer;
    }

    node query(int l, int r) {
        return query(1, 0, n - 1, l, r);
    }

    void update(int index, int idx, int ss, int se) {
        if (idx < ss || idx > se)
            return;
        if (ss == se) {
            single(tree[index], ss);
            return;
        }
        int mid = (ss + se) / 2;
        update(2 * index, idx, ss, mid);
        update(2 * index + 1, idx, mid + 1, se);
        merge(tree[index], tree[2 * index], tree[2 * index + 1]);
    }

    void update(int idx,int delta) {
        a[idx]=delta;
        update(1, idx, 0, n - 1);
    }
};

/*-----------SGT___with_lazy------*/
int a[200001];
//for range min query and range update
class segte {
public:
    struct node{ //req variable
        int val; //default value
        node() {
            val=INT_MAX;
        }
    };
    int n{};
    vector<node> tree;
    vector<int>lazy;
    vector<bool>prsnt;
    node neutral;

    void init(int N) {
        n = N;
        tree.resize(4 * n + 1);
        //default values
        lazy.assign(4*n+1,0);
        prsnt.assign(4*n+1,false);
    }

    //merge function
    static void merge(node &curr, node &left, node &right) {
        curr.val =min( left.val , right.val);
    }


    void prop(int index,int ss,int se)
    {
        tree[index].val+=lazy[index];
        if(ss!=se)
        {
            lazy[2*index]+=lazy[index];
            lazy[2*index+1]+=lazy[index];
            prsnt[2*index]=prsnt[2*index+1]=true;
        }
        lazy[index]=0;
        prsnt[index]=false;
    }

    void build(int index, int ss, int se) {
        if (ss == se) {
            tree[index].val=a[ss];
            return;
        }
        int mid = (ss + se) / 2;
        build(2 * index, ss, mid);
        build(2 * index + 1, mid + 1, se);
        merge(tree[index], tree[2 * index], tree[2 * index + 1]);
    }

    void build() {
        build(1, 0, n - 1);
    }

    node query(int index,int ss,int se,int qs,int qe)
    {
        if(prsnt[index])
            prop(index,ss,se);
        if (qs > se || qe < ss) return neutral;
        if (qs <= ss && qe >= se) return tree[index];
        int mid = (ss + se) / 2;
        node left = query(2 * index, ss, mid, qs, qe);
        node right = query(2 * index + 1, mid + 1, se, qs, qe);
        node mer;
        merge(mer, left, right);
        return mer;
    }
    node query(int l, int r) {
        return query(1, 0, n - 1, l, r);
    }
    void update(int index, int ss, int se,int l,int r,int inc) {
        if (prsnt[index])
            prop(index, ss, se);
        if (r < ss || l > se)
            return;
        if (ss >= l && se <= r) {
            prsnt[index] = true;
            lazy[index] += inc;
            prop(index, ss, se);
            return;
        }
        int mid = (ss + se) / 2;
        update(2 * index, ss, mid, l, r, inc);
        update(2 * index + 1, mid + 1, se, l, r, inc);
        merge(tree[index], tree[2 * index], tree[2 * index + 1]);
    }
    void update(int l,int r,int inc) {
        update(1, 0, n - 1,l,r,inc);
    }
};
/*------------2D_SGT--------------*/
class segte2d{
public:
    int n,m;
    vector<vector<int>>tr;
    vector<vector<int>>a;
    segte2d(int _n,int _m,vector<vector<int>>&_a){
        a=_a;
        n=_n;
        m=_m;
        tr.resize(4*n+1,vector<int>(4*m+1,0));
    }

    int f(int x,int y) {
        return __gcd(x,y);
    }
    void build_y(int vx,int vy,int sx,int ex,int sy,int ey){
        if(sy==ey) {
            if (sx == ex) tr[vx][vy] = a[sx][sy];
            else tr[vx][vy] = f(tr[2 * vx + 1][vy] , tr[2 * vx + 2][vy]);
            return;
        }
        int my=(sy+ey)/2;
        build_y(vx,2*vy+1,sx,ex,sy,my);
        build_y(vx,2*vy+2,sx,ex,my+1,ey);
        tr[vx][vy]=f(tr[vx][2*vy+1],tr[vx][2*vy+2]);
    }
    void build_x(int vx,int sx,int ex){
        if(sx!=ex){
            int mx=(sx+ex)/2;
            build_x(2*vx+1,sx,mx);
            build_x(2*vx+2,mx+1,ex);
        }
        build_y(vx,0,sx,ex,0,m-1);
    }
    void build(){
        build_x(0,0,n-1);
    }
    void update_y(int vx,int vy,int sx,int ex,int sy,int ey,int x,int y,int val){
        if(sy==ey){
            if(sx==ex){
                tr[vx][vy]=f(tr[vx][vy],val);
            }else{
                tr[vx][vy]=f(tr[2*vx+1][vy],tr[2*vx+2][vy]);
            }
            return;
        }
        int mid=(sy+ey)/2;
        if(y<=mid)
            update_y(vx,2*vy+1,sx,ex,sy,mid,x,y,val);
        else
            update_y(vx,2*vy+2,sx,ex,mid+1,ey,x,y,val);
        tr[vx][vy]=f(tr[vx][2*vy+1],tr[vx][2*vy+2]);
    }
    void update_x(int vx,int sx,int ex,int x,int y,int val){
        if(sx!=ex){
            int mid=(sx+ex)/2;
            if(x<=mid){
                update_x(2*vx+1,sx,mid,x,y,val);
            }else{
                update_x(2*vx+2,mid+1,ex,x,y,val);
            }
        }
        update_y(vx,0,sx,ex,0,m-1,x,y,val);
    }
    void update(int x,int y,int val){
        update_x(0,0,n-1,x,y,val);
    }
    int query_y(int vx,int vy,int sy,int ey,int qsy,int qey){
        if(qey<sy || qsy>ey){
            return 0;
        }
        if(sy>=qsy && ey<=qey){
            return tr[vx][vy];
        }
        int my=(sy+ey)/2;
        int A=query_y(vx,2*vy+1,sy,my,qsy,qey);
        int B=query_y(vx,2*vy+2,my+1,ey,qsy,qey);
        return f(A,B);
    }
    int query_x(int vx,int qsx,int qex,int qsy,int qey,int sx,int ex){
        if(qsx>ex || qex<sx)
            return 0;
        if(sx>=qsx && ex<=qex){
            return query_y(vx,0,0,m-1,qsy,qey);
        }
        int mx=(sx+ex)/2;
        int A=query_x(2*vx+1,qsx,qex,qsy,qey,sx,mx) ;
        int B=query_x(2*vx+2,qsx,qex,qsy,qey,mx+1,ex);
        return f(A,B);
    }
    int query(int sx,int ex,int sy,int ey){
        return query_x(0,sx,ex,sy,ey,0,n-1);
    }
};/// my template SGT :
struct item
{
    ll sum,gcd;
    item() :sum((ll)0),gcd((ll)0) {};
    item(ll x) :sum(x),gcd(x) {};
    item(ll x,ll y) :sum(x),gcd(y) {};
};
struct SGT{
    item comb(item&I1,item&I2){
        return item(I1.sum+I2.sum,__gcd(I1.gcd,I2.gcd));
    }
    ll n;
    vector<item> tree;
    vll lazy1,lazy2;
    SGT(ll size) : n(size), tree(4 * size), lazy1(4 * size,1),lazy2(4 * size,1) {}
    void build(const vll &data, ll node, ll start, ll end) {
        if (start == end) {
            tree[node]=item(data[start]);
        } else {
            ll mid = (start + end) / 2;
            build(data, 2 * node, start, mid);
            build(data, 2 * node + 1, mid + 1, end);
            tree[node] =comb(tree[2 * node],tree[2 * node + 1]);
        }
    }
    void push(ll node, ll start, ll end) {
        if (lazy1[node] !=1 || lazy2[node]!=1) {
            ll d=__gcd(lazy1[node],lazy2[node]);
            lazy1[node]/=d;
            lazy2[node]/=d;
            tree[node].sum/=lazy2[node];
            tree[node].gcd/=lazy2[node];
            tree[node].gcd*=lazy1[node];
            tree[node].sum*=lazy1[node];
            if (start != end) {
                lazy1[2 * node]*=lazy1[node];
                lazy1[2 * node + 1]*=lazy1[node];
                lazy2[2 * node]*=lazy2[node];
                lazy2[2 * node + 1]*=lazy2[node];
            }
            lazy1[node]=lazy2[node]=1;
        }
    }
    void updateRange(ll l, ll r, ll value, ll node, ll start, ll end,ll t) {
        push(node, start, end);
        if (start > r || end < l) return;
        if (start >= l && end <= r) {
            if(t==0) lazy1[node]*= value;
            else lazy2[node]*=value;
            push(node, start, end);
            return;
        }
        ll mid = (start + end) / 2;
        updateRange(l, r, value, 2 * node, start, mid,t);
        updateRange(l, r, value, 2 * node + 1, mid + 1, end,t);
        tree[node] =comb(tree[2 * node],tree[2 * node + 1]);
    }
    item query(ll l, ll r, ll node, ll start, ll end) {
        push(node, start, end);
        if (start > r || end < l) return item();
        if (start >= l && end <= r) return tree[node];
        ll mid = (start + end) / 2;
        item it1=query(l, r, 2 * node, start, mid);
        item it2=query(l, r, 2 * node + 1, mid + 1, end);
        return comb(it1,it2);
    }
    void updateElement(ll idx, ll value, ll node, ll start, ll end) {
        push(node, start, end);
        if (start == end) {
            tree[node] =item(value);
        } else {
            ll mid = (start + end) / 2;
            if (idx <= mid) {
                updateElement(idx, value, 2 * node, start, mid);
            } else {
                updateElement(idx, value, 2 * node + 1, mid + 1, end);
            }
            tree[node] =comb(tree[2 * node],tree[2 * node + 1]);
        }
    }
    item queryElement(ll idx, ll node, ll start, ll end) {
        push(node, start, end);
        if (start == end) return tree[node];
        ll mid = (start + end) / 2;
        if (idx <= mid) return queryElement(idx, 2 * node, start, mid);
        else return queryElement(idx, 2 * node + 1, mid + 1, end);
    }
    void updateRange(ll l, ll r, ll value,ll t) {
        updateRange(l, r, value, 1, 0, n - 1,t);
    }
    item query(ll l, ll r) {
        return query(l, r, 1, 0, n - 1);
    }
    void updateElement(ll idx, ll value) {
        updateElement(idx, value, 1, 0, n - 1);
    }
    item queryElement(ll idx) {
        return queryElement(idx, 1, 0, n - 1);
    }
    void build(const vll &data) {
        build(data, 1, 0, n - 1);
    }
};
/*-------------sparse table1D---------*/
class SparseTable {
private:
    static const int N = 200000;
    static const int M = 21;
    int tab[N+1][M+1];
    int L[N+1];
    int n;
    int f(int x, int y) const {
        return max(x, y);
    }
public:
    SparseTable(int _n,vector<int>&_a) : n(_n) {
        for (int i = 2; i <= n; ++i) {
            L[i] = L[i / 2] + 1;
        }
        for (int i = 0; i < n; ++i) {
            tab[i][0] =_a[i];
        }
    }
    void build() {
        for (int j = 1; j <= M; ++j) {
            for (int i = 0; i + (1 << j) - 1 < n; ++i) {
                tab[i][j] = f(tab[i][j-1], tab[i + (1 << (j-1))][j-1]);
            }
        }
    }
    int query(int l, int r) const {
        int len = r - l + 1;
        int idx = l;
        int tot = 0; // Neutral 
        for (int j = M; j >= 0; --j) {
            if (len & (1 << j)) {
                tot = f(tot, tab[idx][j]);
                idx += (1 << j);
            }
        }
        return tot;
    }
    int query_i(int l, int r) const {
        int lg = L[r - l + 1];
        return f(tab[l][lg], tab[r - (1 << lg) + 1][lg]);
    }
};

class SparseTable2D {
private:
    static const int N = 1000;
    static const int M = 10;
    int tab[N+1][N+1][M+1][M+1];
    int L[N+1];
    int n, m;
    int f(int x, int y, int z = 0, int w = 0) const {
        return max({x, y, z, w});
    }
public:
    SparseTable2D(int _n, int _m,vector<vector<int>>&_a) : n(_n), m(_m) {
        for (int i = 2; i <= N; ++i) {
            L[i] = L[i / 2] + 1;
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                tab[i][j][0][0] =_a[i][j];
            }
        }
    }
    void build() {
        for (int x = 1; x <= M; ++x) {
            for (int i = 0; i + (1 << x) - 1 < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    tab[i][j][x][0] = f(tab[i][j][x-1][0], tab[i + (1 << (x-1))][j][x-1][0]);
                }
            }
        }
        for (int y = 1; y <= M; ++y) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j + (1 << y) - 1 < m; ++j) {
                    tab[i][j][0][y] = f(tab[i][j][0][y-1], tab[i][j + (1 << (y-1))][0][y-1]);
                }
            }
        }
        for (int x = 1; x <= M; ++x) {
            for (int y = 1; y <= M; ++y) {
                for (int i = 0; i + (1 << x) - 1 < n; ++i) {
                    for (int j = 0; j + (1 << y) - 1 < m; ++j) {
                        tab[i][j][x][y] = f(
                            tab[i][j][x-1][y-1],
                            tab[i + (1 << (x-1))][j][x-1][y-1],
                            tab[i][j + (1 << (y-1))][x-1][y-1],
                            tab[i + (1 << (x-1))][j + (1 << (y-1))][x-1][y-1]
                        );
                    }
                }
            }
        }
    }
    int query_i(int x, int y, int _x, int _y) const {
        if (_x < x || _y < y) return 0; // Neutral 
        int lx = L[_x - x + 1];
        int ly = L[_y - y + 1];
        return f(
            tab[x][y][lx][ly],
            tab[_x - (1 << lx) + 1][y][lx][ly],
            tab[x][_y - (1 << ly) + 1][lx][ly],
            tab[_x - (1 << lx) + 1][_y - (1 << ly) + 1][lx][ly]
        );
    }
};

/*--------MO--------------*/
struct Mo {
    ll n, q, block_size, current_answer;
    vll arr, freq, answers;
    vll compressed_arr, compressed;
    vector<unordered_map<ll,ll>> compress_map;
    vll decompress;
    struct Query {
        ll l, r, idx;
    };
    vector<Query> queries;
    
    Mo(ll n,ll q) : n(n), q(q),block_size((ll)sqrt(n)),current_answer(0), arr(n), freq(n), answers(q) {
        queries.resize(q);
    }

    // Coordinate compression
    void compressArray() {
        compressed = arr;
        sort(compressed.begin(), compressed.end());
        compressed.erase(unique(compressed.begin(), compressed.end()), compressed.end());
        
        unordered_map<ll,ll> compress_map;
        for (ll i = 0; i < compressed.size(); i++) {
            compress_map[compressed[i]] = i;
        }

        // Convert the original array to compressed array
        compressed_arr.resize(n);
        decompress.resize(n);
        for (ll i = 0; i < n; i++) {
            compressed_arr[i] = compress_map[arr[i]];
            decompress[compressed_arr[i]]=arr[i];
        }
    }

    // Mo's algorithm comparison function
    static bool mo_cmp(const Query &a, const Query &b,ll block_size) {
        if (a.l / block_size != b.l / block_size)
            return a.l / block_size < b.l / block_size;
        return a.r < b.r;
    }

    // Process all queries using Mo's algorithm
    void processQueries() {
        sort(queries.begin(), queries.end(), [this](const Query &a, const Query &b) {
            return mo_cmp(a, b, block_size);
        });

        freq.assign(compressed.size(), 0);
        current_answer = 0;
        ll curr_l = 0, curr_r = -1;

        auto add = [&](ll x) {
            freq[x]++;
            if (freq[x]%2) current_answer++;
            else current_answer--;
        };

        auto remove = [&](int x) {
            if (freq[x]%2) current_answer--;
            else current_answer++;
            freq[x]--;
        };

        // Process each query
        for (const Query &query : queries) {
            ll l = query.l;
            ll r = query.r;

            while (curr_r < r) {
                curr_r++;
                add(compressed_arr[curr_r]);
            }
            while (curr_r > r) {
                remove(compressed_arr[curr_r]);
                curr_r--;
            }
            while (curr_l < l) {
                remove(compressed_arr[curr_l]);
                curr_l++;
            }
            while (curr_l > l) {
                curr_l--;
                add(compressed_arr[curr_l]);
            }

            answers[query.idx] =(current_answer==0);
        }
    }

    // Add a query
    void addQuery(ll l,ll r,ll idx) {
        queries[idx] = {l - 1, r - 1, idx}; // Assuming 0-based indexing
    }

    // Input function
    void inputArray() {
        for (ll i = 0; i < n; i++) {
            cin >> arr[i];
        }
    }

    void inputQueries() {
        for (ll i = 0; i < q; i++) {
            ll l, r;
            cin >> l >> r;
            addQuery(l, r, i);
        }
    }

    // Output results
    void outputAnswers() {
        for (ll i = 0; i < q; i++) {
            if(answers[i]) YES
            else NO
        }
    }
};

void sol() {
    ll n, m;
    cin >> n >> m;
    Mo mo(n,m);
    mo.inputArray();
    mo.inputQueries();
    mo.compressArray();
    mo.processQueries();
    mo.outputAnswers();
}


/*-------DSU--------------*/
struct DSU {
    vll sz, parent;
    ll Nbrset;
    // Initialize DSU with 'n' elements
    void Init(ll n) {
        // Resize the vectors to hold 'n' elements
        sz.resize(n);
        parent.resize(n);
        // Initialize each element to be its own parent (each element is its own set)
        for (ll i = 0; i < n; i++) {
            parent[i] = i;
            sz[i] = 1; // Each set initially has a size of 1
        }
        Nbrset=n;
    }

    // Find the root of the set containing 'u' with path compression
    ll get(ll u) {
        // If 'u' is not its own parent, recursively find the root and compress the path
        return (parent[u] == u) ? u : (parent[u] = get(parent[u]));
    }

    // Union the sets containing 'u' and 'v'
    bool unite(ll u,ll v) {
        // Find the roots of the sets containing 'u' and 'v'
        u = get(u);
        v = get(v);
        // If they are already in the same set, return false
        if (u == v) return false;
        // Union by size: attach the smaller tree under the root of the larger tree
        if (sz[u] < sz[v]) swap(u, v);
        parent[v] = u; // Make 'u' the parent of 'v'
        sz[u] += sz[v]; // Update the size of the root
        Nbrset--;
        return true; // Indicate that a union was performed
    }
    ll  numDisjointSets() { return Nbrset; }
    // Check if 'u' and 'v' are in the same set
    bool same_set(ll u,ll v) {
        return get(u) == get(v);
    }

    // Get the size of the set containing 'u'
    ll size(ll u) {
        return sz[get(u)];
    }
};

class dsu { /// with_role_back
public:
    vector<int> parent;
    vector<int> size;
    stack<tuple<int, int, int, int>> stk;
    stack<int> flag;
    int comp = 0;

    explicit dsu(int a) {
        parent.resize(a);
        size.resize(a);
        comp = a;
        for (int i = 0; i < a; i++) {
            parent[i] = i;
            size[i] = 1;
        }
    }

    int par(int i) {
        if (i == parent[i])
            return i;
        return par(parent[i]);
    }

    void unite(int a, int b) {
        a = par(a);
        b = par(b);
        if (a != b) {
            comp--;
            stk.push({a, b, size[a], parent[a]});
            parent[a] = b;
            size[b] += size[a];
            size[a] = 0;
        }
    }

    void save() {
        flag.push(stk.size());
    }

    void roll_back() {
        int a = get<0>(stk.top());
        int b = get<1>(stk.top());
        int sz = get<2>(stk.top());
        int rp = get<3>(stk.top());
        stk.pop();
        comp++;
        parent[a] = rp;
        size[a] = sz;
        size[b] -= size[a];
    }

    void roll() {
        int last = flag.top();
        flag.pop();
        while (stk.size() > last)
            roll_back();
    }
};

//1-base
class BIT {
public:
    int n;
    vector<int> tree;

    void init(int _n) {
        n = _n;
        tree.resize(n + 5, 0);
    }

    void update(int idx, int inc) {
        for (int i = idx; i <= n; i += i & -i)
            tree[i] += inc;
    }

    int query(int idx) {
        int sum = 0;
        assert(idx<=n);
        for (int i = idx; i >= 1; i -= i & -i)
            sum += tree[i];
        return sum;
    }
};

class BIT_update_range_query {
public:
    int n;
    vector<int> tree[2];
    void init(int _n) {
        n = _n;
        tree[0].resize(n + 5, 0);
        tree[1].resize(n + 5, 0);
    }
    void upd(int which, int idx, int inc) {
        for (int j = idx; j <= n; j += j & -j) {
            tree[which][j] += inc;
        }
    }
    void range_upd(int l, int r, int x) {
        upd(0, l, x);
        upd(0, r + 1, -x);
        upd(1, l, x * (l - 1));
        upd(1, r + 1, -x * r);
    }
    int qry(int which, int idx) {
        int tot = 0;
        for (int j = idx; j <= n; j += j & -j) {
            tot += tree[which][j];
        }
        return tot;
    }
    int prefix_qry(int idx) {
        return qry(0, idx) * idx - qry(1, idx);
    }
};


struct sd{
    vector<int>a,bv,bno,bl,br;
    int n,bsz,bcnt;

    sd(int _n,vector<int>&_a){
        n=_n;
        bsz=sqrt(n);
        bcnt=(n+bsz-1)/bsz;

        a=_a;
        bno.resize(n,0);

        bl.resize(bcnt,-1);
        br.resize(bcnt,-1);
        bv.resize(bcnt,0);

        for(int i =0;i<n;i++){
            bno[i]=i/bsz;
            if(bl[bno[i]]==-1) bl[bno[i]]=i;
            br[bno[i]]=i;

            bv[bno[i]]+=a[i];
        }
    }

    void upd(int i,int x){
        bv[bno[i]]-=a[i];
        a[i]=x;
        bv[bno[i]]+=a[i];
    }

    int qry(int l,int r){
        int tot=0;
        if(bno[l]==bno[r]){
            for(int j=l;j<=r;j++){
                tot+=a[j];
            }
        }else{
            for(int j=l;j<=br[bno[l]];j++){
                tot+=a[j];
            }
            for(int j=bno[l]+1;j<bno[r];j++){
                tot+=bv[j];
            }
            for(int j=bl[bno[r]];j<=r;j++){
                tot+=a[j];
            }
        }
        return tot;
    }

};


struct sd_lazy{
    vector<int>a,bv,bno,bl,br,lazy;
    int n,bsz,bcnt;

    sd_lazy(int _n,vector<int>&_a){
        n=_n;
        bsz=sqrt(n);
        bcnt=(n+bsz-1)/bsz;

        a=_a;
        bno.resize(n,0);

        bl.resize(bcnt,-1);
        br.resize(bcnt,-1);
        bv.resize(bcnt,0);
        lazy.resize(bcnt,0);

        for(int i =0;i<n;i++){
            bno[i]=i/bsz;
            if(bl[bno[i]]==-1) bl[bno[i]]=i;
            br[bno[i]]=i;

            bv[bno[i]]+=a[i];
        }
    }

    void prop(int block_num){
        
        if(lazy[block_num]==0){
            return;
        }
        for(int j=bl[block_num];j<=br[block_num];j++){
            a[j]+=lazy[block_num];
        }
        bv[block_num]+=lazy[block_num]*(br[block_num]-bl[block_num]+1);
        lazy[block_num]=0;
    }

    void upd(int l,int r,int x){
        if(bno[l]==bno[r]){
            prop(bno[l]);
            for(int j=l;j<=r;j++){
                a[j]+=x;
                bv[bno[l]]+=x;
            }
        }else{
            prop(bno[l]);
            for(int j=l;j<=br[bno[l]];j++){
                a[j]+=x;
                bv[bno[l]]+=x;
            }
            for(int j=bno[l]+1;j<bno[r];j++){
                lazy[j]+=x;
            }
            prop(bno[r]);
            for(int j=bl[bno[r]];j<=r;j++){
                a[j]+=x;
                bv[bno[r]]+=x;
            }
        }
    }

    int qry(int l,int r){
        int tot=0;
        if(bno[l]==bno[r]){
            prop(bno[l]);
            for(int j=l;j<=r;j++){
                tot+=a[j];
            }
        }else{
            prop(bno[l]);
            for(int j=l;j<=br[bno[l]];j++){
                tot+=a[j];
            }
            for(int j=bno[l]+1;j<bno[r];j++){
                tot+=bv[j];
                tot+=lazy[j]*(br[j]-bl[j]+1);
            }
            prop(bno[r]);
            for(int j=bl[bno[r]];j<=r;j++){
                tot+=a[j];
            }
        }
        return tot;
    }

};


/// Treap
const int N = 2e5 + 9;
mt19937 rnd(chrono::steady_clock::now().time_since_epoch().count());
struct node {
  node *l, *r;
  int key, prior;
  node(int id) {
    l = r = nullptr;
    key = id;
    prior = rnd();
  }
};
struct treap {
  node *root;
  treap() {
    root = nullptr;
  }
  void split(node *t, int pos, node *&l, node *&r) {
    if (t == nullptr) {
      l = r = nullptr;
      return;
    }
    if (t->key <= pos) {
      split(t->r, pos, l, r);
      t->r = l;
      l = t;
    } else {
      split(t->l, pos, l, r);
      t->l = r;
      r = t;
    }
  }
  node* merge(node *l, node *r) {
    if (!l || !r) return l ? l : r;
    if (l->prior < r->prior) {
      l->r = merge(l->r, r);
      return l;
    }
    r->l = merge(l, r->l);
    return r;
  }
  node* merge_treap(node *l, node *r) {
    if(!l) return r;
    if(!r) return l;
    if(l->prior < r->prior) swap(l, r);
    node *L, *R;
    split(r, l->key, L, R);
    l->r = merge_treap(l->r, R);
    l->l = merge_treap(L, l->l);
    return l;
  }
  void insert(int id) {
    node *l, *r;
    split(root, id, l, r);
    root = merge(merge(l, new node(id)), r);
  }
  node* erase(int L, int R) {
    node *l, *r, *mid, *mr;
    split(root, L - 1, l, r);
    split(r, R, mid, mr);
    root = merge(l, mr);
    return mid;
  }
  void combine(node *x) { //combine with another treap
    root = merge_treap(root, x);
  }
  vector<int> ans;
  void dfs(node* cur) {
    if(!cur) return;
    ans.push_back(cur -> key);
    dfs(cur -> l);
    dfs(cur -> r);
  }
  vector<int> get() {
    ans.clear();
    dfs(root);
    return ans;
  }
} t[N];
int ans[N];
int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  int n;
  cin >> n;
  for(int i = 1; i <= n; i++) {
    int k;
    cin >> k;
    t[k].insert(i);
  }
  int q;
  cin >> q;
  while(q--) {
    int l, r, x, y;
    cin >> l >> r >> x >> y;
    node *p = t[x].erase(l, r);
    t[y].combine(p);
  }
  for(int i = 1; i < N; i++) {
    auto p = t[i].get();
    for(auto x : p) ans[x] = i;
  }
  for(int i = 1; i <= n; i++) cout << ans[i] << " \n"[i == n];
  return 0;
}
// https://codeforces.com/problemset/problem/911/G


/// Trie
struct Trie {
  static const int B = 31;
  struct node {
    node* nxt[2];
    int sz;
    node() {
      nxt[0] = nxt[1] = NULL;
      sz = 0;
    }
  }*root;
  Trie() {
    root = new node();
  }
  void insert(int val) {
    node* cur = root;
    cur -> sz++;
    for (int i = B - 1; i >= 0; i--) {
      int b = val >> i & 1;
      if (cur -> nxt[b] == NULL) cur -> nxt[b] = new node();
      cur = cur -> nxt[b];
      cur -> sz++;
    }
  }
  int query(int x, int k) { // number of values s.t. val ^ x < k
    node* cur = root;
    int ans = 0;
    for (int i = B - 1; i >= 0; i--) {
      if (cur == NULL) break;
      int b1 = x >> i & 1, b2 = k >> i & 1;
      if (b2 == 1) {
        if (cur -> nxt[b1]) ans += cur -> nxt[b1] -> sz;
        cur = cur -> nxt[!b1];
      } else cur = cur -> nxt[b1];
    }
    return ans;
  }
  int get_max(int x) { // returns maximum of val ^ x
    node* cur = root;
    int ans = 0;
    for (int i = B - 1; i >= 0; i--) {
      int k = x >> i & 1;
      if (cur -> nxt[!k]) cur = cur -> nxt[!k], ans <<= 1, ans++;
      else cur = cur -> nxt[k], ans <<= 1;
    }
    return ans;
  }
  int get_min(int x) { // returns minimum of val ^ x
    node* cur = root;
    int ans = 0;
    for (int i = B - 1; i >= 0; i--) {
      int k = x >> i & 1;
      if (cur -> nxt[k]) cur = cur -> nxt[k], ans <<= 1;
      else cur = cur -> nxt[!k], ans <<= 1, ans++;
    }
    return ans;
  }
  void del(node* cur) {
    for (int i = 0; i < 2; i++) if (cur -> nxt[i]) del(cur -> nxt[i]);
    delete(cur);
  }
} t;
// int16_t main() {
//   ios_base::sync_with_stdio(0);
//   cin.tie(0);
//   int n, k;
//   cin >> n >> k;
//   int cur = 0;
//   long long ans = 1LL * n * (n + 1) / 2;
//   t.insert(cur);
//   for (int i = 0; i < n; i++) {
//     int x;
//     cin >> x;
//     cur ^= x;
//     ans -= t.query(cur, k);
//     t.insert(cur);
//   }
//   cout << ans << '\n';
//   return 0;
// }

