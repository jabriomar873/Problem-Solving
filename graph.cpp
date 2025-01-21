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

struct Edge {
    ll u, v, weight;
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
  
};

/////////--------------Tree-----------------//////

///______LCA:


const ll LOG =18; // log2(MAXN)+1
vll adj[MAXN]; 
ll up[MAXN][LOG]; // up[v][j] is the 2^j-th Ac of node v
ll depth[MAXN]; 
ll sz[MAXN];
ll A[MAXN];
ll st[MAXN];
ll en[MAXN];
ll timer = 0;
void dfsLCA(ll v,ll p) {
    up[v][0] = p;
    for (ll i = 1; i < LOG; ++i) {
        if (up[v][i-1] != -1) {
            up[v][i] = up[up[v][i-1]][i-1];
        } else {
            up[v][i] = -1;
        }
    }
    for (ll u : adj[v]) {
        if (u != p) {
            depth[u] = depth[v] + 1;
            dfsLCA(u, v);
        }
    }
}
void preprocess(ll root=1) {
    ////memset(up,-1,sizeof(up));
    depth[root]=0;
    dfsLCA(root,-1);////-1 is the father of root 
}
ll find_kth_ancestor(ll x,ll k) {
    for (ll j = 0; j < LOG; ++j) {
        if ((k >> j) & 1) {
            x = up[x][j];
            if (x == -1) break;
        }
    }
    return x;
}
ll LCA(ll u,ll v) {
    if (depth[u] < depth[v]) swap(u, v);
    ll diff = depth[u] - depth[v];
    for (ll i = 0; i < LOG; ++i) {
        if ((diff >> i) & 1) {
            u = up[u][i];
        }
    }
    if (u == v) return u;
    for (ll i = LOG - 1; i >= 0; --i) {
        if (up[u][i] != up[v][i]) {
            u = up[u][i];
            v = up[v][i];
        }
    }
    return up[u][0];
}
void euler_tour(ll x, ll p) {
    st[x] =timer;
    timer++;
    ///cout<<x<<"----st="<<st[x]<<endl;
    for (const ll &e : adj[x])
        if (e != p) euler_tour(e,x);
    en[x] =timer-1;
    ///cout<<x<<"---en="<<en[x]<<endl;
}




void findCycleUtil(ll v,const vector<vll>& adj,vll& visited,vll& parent,vbo& inPath,vll& cycle,bool &cycleFound) {
    visited[v] = 1; // mark as visited
    inPath[v] = true; // mark as in the current path
    
    for (ll u:adj[v]) {
        if (!visited[u]) { // if not visited
            parent[u] = v;
            findCycleUtil(u,adj,visited,parent,inPath,cycle,cycleFound);
            if (cycleFound) return; // if cycle is found, return immediately
        } else if (inPath[u]) { // if in current path, cycle detected
            cycleFound=true;
            // Trace back to get the cycle
            for (ll x =v; x !=u; x=parent[x]) {
                cycle.push_back(x);
            }
            cycle.push_back(u);
            cycle.push_back(v);
            return;
        }
    }
    inPath[v]=false; // unmark from the current path
}

vector<ll> findCycle(ll N,const vector<vll>& adj) {
    vll visited(N + 1,0);
    vll parent(N + 1,-1);
    vector<bool> inPath(N + 1,false);
    vll cycle;
    bool cycleFound = false;
    for (ll v=1;v<= N;v++) {
        if (!visited[v]) {
            findCycleUtil(v,adj,visited,parent,inPath,cycle,cycleFound);
            if (cycleFound) {
                reverse(cycle.begin(), cycle.end());
                break;
            }
        }
    }
    return cycle; /// if cycle empty =there is no cycle
}   ///// cycle= 1 3 5 6 1 comme ça

void dffs(ll v,vector<bool>& visited,stack<ll>& Stack,vector<vll>& adj) {
    visited[v] = true;

    for (int i : adj[v]) {
        if (!visited[i]) {
            dffs(i, visited, Stack, adj);
        }
    }

    Stack.push(v);
}

vector<ll> topologicalSort(ll V, vector<vll>& adj) {
    stack<ll> Stack;
    vector<bool> visited(V, false);
        //!!!!!! V is not the rote or any thing just the size          
                                       //// using DFS topological sort
    for (int i = 0; i < V; i++) {
        if (!visited[i]) {
            dffs(i, visited,Stack, adj);
        }
    }

    vll topo_order;
    while (!Stack.empty()) {
        topo_order.push_back(Stack.top());
        Stack.pop();
    }

    return topo_order;
}
void find_cycles(const vll& p, vector<vll>& cycles) {   /// for permuation
    ll n = p.size();
    vector<bool> visited(n, false);
    for (ll i = 0; i < n; ++i) {
        if (!visited[i]) {
            vll cycle;
            ll x = i;
            while (!visited[x]) {
                visited[x] = true;
                cycle.push_back(x);
                x = p[x] - 1;
            }
            cycles.push_back(cycle);
        }
    }
}


pair<vector<ll>, vector<ll>> dijkstra(const vector<vector<pll>>& graph, ll source) {
    ll n = graph.size(); // Number of vertices
    vector<ll> dist(n, INF); // Initialize distances to infinity
    vector<ll> pred(n, -1); // To store the predecessor of each vertex
    priority_queue<pll, vector<pll>, greater<pll>> pq; // Min-heap for priority queue

    dist[source] = 0; // Distance from source to itself is 0
    pq.push({0, source}); // Push the source vertex into the priority queue

    while (!pq.empty()) {
        ll u = pq.top().second; // Extract vertex with minimum distance
        ll d = pq.top().first; // Extracted distance
        pq.pop();

        // If the extracted distance is greater than the currently known distance, skip
        if (d > dist[u]) continue;

        // Check all neighbors of u
        for (const auto& neighbor : graph[u]) {
            ll v = neighbor.first; // Neighbor vertex
            ll w = neighbor.second; // Edge weight from u to v

            // Relaxation step
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w; // Update distance if a shorter path is found
                pred[v] = u; // Set the predecessor of v to u
                pq.push({dist[v], v}); // Push the updated distance and vertex to the priority queue
            }
        }
    }

    return {dist, pred}; // Return both the distance and predecessor vectors
}

vector<ll> get_shortest_path(ll source, ll target, const vector<ll>& pred) {
    vector<ll> path;
    for (ll at = target; at != -1; at = pred[at]) {
        path.push_back(at);
    }
    reverse(path.begin(), path.end());
    if (path[0] != source) {
        return {}; // Return an empty path if there is no path from source to target
    }
    return path;
}


struct edge
{
   ll I,E,W;
   edge():I(0LL),E(0LL),W(0LL){};
   edge(ll u,ll v,ll w):I(u),E(v),W(w){};
};
pair<bool,vector<ll>> bellmanFord(vector<edge> graph,ll vertex,ll Sr,ll Ds ){
    vector<ll> distances(vertex+1,-INF);
    distances[Sr] = 0;
    for(ll i=0;i<vertex-1;i++){
        for(ll j=0;j<graph.size();j++){
            if(distances[graph[j].E] <distances[graph[j].I] + graph[j].W){ ///< pour max > pour min 
                distances[graph[j].E] = distances[graph[j].I] + graph[j].W;
            }
        }
    }
    for(ll j=0;j<graph.size();j++){
        if(distances[graph[j].E] <distances[graph[j].I] + graph[j].W){ ///< pour max > pour min 
            return make_pair(false,vector<ll>());
        }
    }
    return make_pair(true,distances);
}

/*-----bridges+point_articaltion-----*/
const int N=200001;
vector<int>gr[N];
int dis[N];
int low[N];
int tme;
set<int>art_p;
vector<pair<int,int>>bridge;
void dfs(int node,int par) {
    dis[node] = low[node] = tme++;
    int child = 0;
    for (auto c:gr[node]) {
        if (!dis[c]) {
            child++;
            dfs(c, node);
            low[node] = min(low[node], low[c]);
            if (par != -1 && low[c] >= dis[node]) {
                art_p.insert(node);
            }
            if (low[c] > dis[node]) {
                bridge.push_back({node, c});
            }
        } else if (c != par) {
            low[node] = min(low[node], dis[c]);
        }
    }
    if (par == -1 && child > 1)
        art_p.insert(1);
}
void solve() {
    int n, m;
    cin >> n >> m;
    for (int i = 0; i < m; i++) {
        int x, y;
        cin >> x >> y;
        --x;
        --y;
        gr[x].push_back(y);
        gr[y].push_back(x);
    }
    tme = 1;
    dfs(0, -1);
    cout << "cut vertices" << endl;
    for (auto c:art_p) {
        cout << c << " ";
    }
    cout << endl;
    cout << "bridges" << endl;
    for (auto c:bridge) {
        cout << c.first << " " << c.second << endl;
    }
}

/*____________graph____________*/
class Graph { //// Briges
public:
    Graph(ll n) : n(n), adj(n), timer(0), bridges_count(0) {}
    void addEdge(ll u,ll v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    ll countBridges() {
        tin.assign(n, -1);
        low.assign(n, -1);
        visited.assign(n, false);
        for (ll i = 0; i < n; i++) {
            if (!visited[i]) {
                dfs(i, -1);
            }
        }
        return bridges_count;
    }
private:
    ll n;
    vector<vll> adj;
    vll tin, low;
    vbo visited;
    ll timer, bridges_count;
    void dfs(ll v,ll parent) {
        visited[v] = true;
        tin[v] = low[v] = timer++;
        for (ll to : adj[v]) {
            if (to == parent) continue;
            if (visited[to]) {
                // Back edge
                low[v] = min(low[v], tin[to]);
            } else {
                // Tree edge
                dfs(to, v);
                low[v] = min(low[v], low[to]);

                if (low[to] > tin[v]) {
                    bridges_count++;
                }
            }
        }
    }
};


/// Floyed
vector<vll> floydWarshall(ll n,vector<pair<pll,ll>>& edges) { /// {{a-b}-c} adge from a to b weighted c
    vector<vll> dist(n,vll(n,INF));             /// !!! if the graph is directed 
                                                ///   so dist[a][b]!=dist[b][a]
    // Initialize diagonal elements to 0
    for (ll i = 0; i < n; ++i) {
        dist[i][i] = 0;
    }
    // Update distances for existing edges
    for (auto& edge : edges) {
        ll a = edge.first.first-1;
        ll b = edge.first.second-1;
        long long c = edge.second;
        dist[a][b] = min(dist[a][b],c);
        dist[b][a] = min(dist[b][a],c);
    }
    // Floyd-Warshall algorithm
    for (ll k = 0; k < n; ++k) {
        for (ll i = 0; i < n; ++i) {
            for (ll j =0; j < n; ++j) {
                if (dist[i][k] < INF && dist[k][j] < INF) {
                    dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
                }
            }
        }
    }
    return dist;
}


/*------strongly connected component SCC-------------*/

vector<vll> adj, adjT;
bool visited[MAXN];
vll component[MAXN];
ll nbrconnected;
stack<ll> stac;
mpll mp;

void dfs1(ll node) {
    visited[node] = true;
    for (auto v : adj[node]) {
        if (!visited[v]) {
            dfs1(v);
        }
    }
    stac.push(node);
}

void dfs2(ll node) {
    visited[node] = true;
    component[nbrconnected].push_back(node);
    mp[node] = nbrconnected;
    for (auto v : adjT[node]) {
        if (!visited[v]) {
            dfs2(v);
        }
    }
}

void findSCCs(ll n) {
    memset(visited, false, sizeof(visited));
    for (ll i = 1; i <= n; ++i) {
        if (!visited[i]) {
            dfs1(i);
        }
    }
    memset(visited, false, sizeof(visited));
    while (!stac.empty()) {
        ll node = stac.top();
        stac.pop();
        if (!visited[node]) {
            nbrconnected++;
            dfs2(node);
        }
    }
}

void sol() {
    ll n;
    cin >> n;
    ll k;
    cin >> k;
    nbrconnected = 0;
    vll p(n), q(n);
    lire(p);
    lire(q);

    ///adj.assign(n + 1, vll());
    ///adjT.assign(n + 1, vll());

    for (ll i = 0; i < n - 1; i++) {
        ll u1 = p[i];
        ll v1 = p[i + 1];
        ll u2 = q[i];
        ll v2 = q[i + 1];
        adj[u1].push_back(v1);
        adj[u2].push_back(v2);
        adjT[v1].push_back(u1);
        adjT[v2].push_back(u2);
    }

    findSCCs(n);
    
    if (nbrconnected < k) {
        cout<<"NO";
        return;
    }
    cout<<"Yes"<<endl;
    queue<ll> qiqo;
    sll soso;
    for (ll i = 1; i <= n; i++) {
        ll nbo = mp[p[i - 1]];
        if (soso.count(nbo) == 0) {
            qiqo.push(nbo);
            soso.insert(nbo);
        }
    }

    string s(n, 'a');
    char c = 'a';
    while (!qiqo.empty() && c < 'a' + k) {
        ll xx = qiqo.front();
        qiqo.pop();
        for (auto u : component[xx]) {
            s[u - 1] = c;
        }
        if(c<'z') c++;
    }
    
    cout << s << endl;
}


/*------LCA in O(1) like sp*/
const int N=200005;
const int M=21;
int tin[N],tout[N],lvl[N],tab[2*N][M],L[2*N];
int tme=-1;
vector<int>g[N];
vector<int>ord;

void eut_3(int node,int p){
    tin[node]=++tme;
    ord.push_back(node);
    for(auto &c:g[node]){
        if(c!=p){
            lvl[c]=lvl[node]+1;
            eut_3(c,node);
            ++tme;
            ord.push_back(node);
        }
    }
    tout[node]=tme;
}
void fill(){
    eut_3(0,-1);
    for(int i=2;i<=ord.size();i++){
        L[i]=L[i/2]+1;
    }
    for(int i=0;i<ord.size();i++){
        tab[i][0]=ord[i];
    }
    for(int j=1;j<M;j++){
        for(int i=0;i+(1<<j)-1<ord.size();i++){
            int x=tab[i][j-1];
            int y=tab[i+(1<<(j-1))][j-1];
            if(lvl[x]<lvl[y]){
                tab[i][j]=x;
            }else{
                tab[i][j]=y;
            }
        }
    }
}
int lca(int u,int v){
    if(tin[u]>tin[v]) swap(u,v);
    int l=tin[u],r=tin[v];
    int lg=L[r-l+1];
    int x=tab[l][lg];
    int y=tab[r-(1<<lg)+1][lg];
    if(lvl[x]<lvl[y]){
        return x;
    }else{
        return y;
    }

}
void solve() {
    int n,q;
    cin>>n>>q;
    for(int i=1;i<n;i++){
        int val;
        cin>>val;
        --val;
        g[val].push_back(i);
        g[i].push_back(val);
    }
    fill();
    while(q--){
        int x,y;
        cin>>x>>y;
        --x;--y;
        cout<<lca(x,y)+1<<endl;
    }
}


//// Hameltonien

bool isSafe(ll v, const vector<vector<ll>>& graph, vector<ll>& path, ll pos) {
    if (graph[path[pos - 1]][v] == 0)
        return false;

    for (ll i = 0; i < pos; ++i)
        if (path[i] == v)
            return false;

    return true;
}

bool hamiltonianUtil(const vector<vector<ll>>& graph, vector<ll>& path, ll pos) {
    if (pos == graph.size()) {
        // Hamiltonian cycle found if graph is complete
        if (graph[path[pos - 1]][path[0]] == 1)
            return true;
        else
            return false;
    }

    for (ll v = 1; v < graph.size(); ++v) {
        if (isSafe(v, graph, path, pos)) {
            path[pos] = v;

            if (hamiltonianUtil(graph, path, pos + 1))
                return true;

            path[pos] = -1;
        }
    }

    return false;
}

bool hamiltonianCycle(const vector<vector<ll>>& graph, vector<ll>& path) {
    path[0] = 0;
    if (!hamiltonianUtil(graph, path, 1))
        return false;

    return true;
}

/// Euler
void dfsE(const vector<vector<ll>>& adj, ll u, vector<bool>& visited, vector<ll>& tempCircuit) {
    visited[u] = true;
    tempCircuit.push_back(u);
    for (ll v : adj[u]) {
        if (!visited[v]) {
            dfsE(adj, v, visited, tempCircuit);
        }
    }
}

void removeEdge(ll u, ll v, vector<vector<ll>>& adj) {
    vector<ll>& adjList = adj[u];
    auto it = find(adjList.begin(), adjList.end(), v);
    if (it != adjList.end()) {
        adjList.erase(it);
    }
}

void eulerAlgorithm(const vector<vector<ll>>& adj, ll u, vector<ll>& circuit) {
    vector<ll> tempCircuit;
    vector<bool> visited(adj.size(), false);
    dfsE(adj, u, visited, tempCircuit);
    for (ll v : tempCircuit) {
        for (ll i = 0; i < adj[v].size(); ++i) {
            ll next = adj[v][i];
            if (next != -1) {
                removeEdge(v, next, const_cast<vector<vector<ll>>&>(adj));  // Remove the visited edge
                circuit.push_back(v);
                break;
            }
        }
    }
    circuit.push_back(u);
}

class Flow {
public:
    ll V;
    vector<vector<ll>> residualGraph;

    Flow(ll V) : V(V) {
        residualGraph.assign(V, vector<ll>(V, 0));
    }

    void addEdge(ll u, ll v, ll capacity) {
        residualGraph[u][v] = capacity;
    }

    bool bfs(ll s, ll t, vector<ll>& parent) {
        vector<bool> visited(V, false);
        queue<ll> q;
        q.push(s);
        visited[s] = true;
        parent[s] = -1;

        while (!q.empty()) {
            ll u = q.front();
            q.pop();

            for (ll v = 0; v < V; v++) {
                if (!visited[v] && residualGraph[u][v] > 0) {
                    q.push(v);
                    parent[v] = u;
                    visited[v] = true;
                }
            }
        }

        return visited[t];
    }

    ll fordFulkerson(ll s, ll t) {
        vector<ll> parent(V, -1);
        ll maxFlow = 0;

        while (bfs(s, t, parent)) {
            ll pathFlow = LLONG_MAX;
            for (ll v = t; v != s; v = parent[v]) {
                ll u = parent[v];
                pathFlow = min(pathFlow, residualGraph[u][v]);
            }

            for (ll v = t; v != s; v = parent[v]) {
                ll u = parent[v];
                residualGraph[u][v] -= pathFlow;
                residualGraph[v][u] += pathFlow;
            }

            maxFlow += pathFlow;
        }

        return maxFlow;
    }

    vpll minCut(ll s, ll t) {
        vector<bool> visited(V, false);
        vector<ll> parent(V, -1);
        queue<ll> q;
        q.push(s);
        visited[s] = true;

        while (!q.empty()) {
            ll u = q.front();
            q.pop();

            for (ll v = 0; v < V; v++) {
                if (!visited[v] && residualGraph[u][v] > 0) {
                    visited[v] = true;
                    q.push(v);
                    parent[v] = u;
                }
            }
        }

        vpll minCutEdges;
        for (ll u = 0; u < V; u++) {
            for (ll v = 0; v < V; v++) {
                if (visited[u] && !visited[v] && residualGraph[u][v] > 0) {
                    minCutEdges.push_back({u, v});
                }
            }
        }

        return minCutEdges;
    }
    void printGraph() {
        cout << "Residual Graph:" << endl;
        for (ll u = 0; u < V; u++) {
            for (ll v = 0; v < V; v++) {
                if (residualGraph[u][v] > 0) {
                    cout << u << " -> " << v << ", Capacity: " << residualGraph[u][v] << endl;
                }
            }
        }
        cout << endl;
    }
};

class Couplage {
public:
    ll V, U;
    vector<vbo> bpGraph;

    Couplage(ll V, ll U) : V(V), U(U) {
        bpGraph.assign(V, vbo(U, false));
    }

    void addEdge(ll u, ll v) {
        bpGraph[u][v] = true;
    }

    bool bpm(ll u, vbo& seen, vll& matchR) {
        for (ll v = 0; v < U; v++) {
            if (bpGraph[u][v] && !seen[v]) {
                seen[v] = true;
                if (matchR[v] < 0 || bpm(matchR[v], seen, matchR)) {
                    matchR[v] = u;
                    return true;
                }
            }
        }
        return false;
    }

    ll maxBPM() {
        vll matchR(U, -1);
        ll result = 0;
        for (ll u = 0; u < V; u++) {
            vbo seen(U, false);
            if (bpm(u, seen, matchR)) {
                result++;
            }
        }
        return result;
    }

    bool isBipartiteUtil(ll src, vll &color) {
        queue<ll> q;
        q.push(src);
        color[src] = 1;

        while (!q.empty()) {
            ll u = q.front();
            q.pop();

            for (ll v = 0; v < U; ++v) {
                if (bpGraph[u][v] && color[v] == -1) {
                    color[v] = 1 - color[u];
                    q.push(v);
                } else if (bpGraph[u][v] && color[v] == color[u]) {
                    return false;
                }
            }
        }
        return true;
    }

    bool isBipartite() {
        vll color(V, -1);
        for (ll i = 0; i < V; i++) {
            if (color[i] == -1) {
                if (!isBipartiteUtil(i, color)) {
                    return false;
                }
            }
        }
        return true;
    }

    vpll printMatching() {
        vll matchR(U, -1);
        vpll matchingPairs;
        for (ll u = 0; u < V; u++) {
            vbo seen(U, false);
            bpm(u, seen, matchR);
        }

        for (ll v = 0; v < U; v++) {
            if (matchR[v] != -1) {
                matchingPairs.push_back({matchR[v], v});
            }
        }
        return matchingPairs;
    }

    vpll findMaximumMatching() {
        vll matchR(U, -1);
        vpll matchingPairs;
        for (ll u = 0; u < V; u++) {
            vbo seen(U, false);
            if (bpm(u, seen, matchR)) {
                for (ll v = 0; v < U; v++) {
                    if (matchR[v] == u) {
                        matchingPairs.push_back({u, v});
                    }
                }
            }
        }
        return matchingPairs;
    }
    
    void addNode(ll u, ll v) {
        bpGraph[u][v] = true;
    }

    bool removeNode(ll u, ll v) {
        if (bpGraph[u][v]) {
            bpGraph[u][v] = false;
            return true;
       
        }
        return false;
    }

    void printGraph() {
        cout << "Bipartite Graph:" << endl;
        for (ll u = 0; u < V; u++) {
            for (ll v = 0; v < U; v++) {
                if (bpGraph[u][v]) {
                    cout << u << " -> " << v << endl;
                }
            }
        }
        cout << endl;
    }
};

