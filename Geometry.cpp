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
 #define read(v) for (auto &u : v) cin >> u;
 #define print(v) for (auto u : v) cout << u << " "; cout <<endl;
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



/*__________Point___________________*/



struct Point {
    ll x, y;
    Point() : x(0), y(0) {};
    Point(ll a, ll b) : x(a), y(b) {};
    friend istream& operator>>(istream& in, Point& p) {
        in >> p.x >> p.y;
        return in;
    }
    friend ostream& operator<<(ostream& out, const Point& p) {
        out << p.x << " " << p.y;
        return out;
    }
    Point operator+(const Point& q) { return Point(x + q.x, y + q.y); }
    Point operator-(const Point& q) { return Point(x - q.x, y - q.y); }
    void operator++() { x++; y++; }
    void operator--() { x--; y--; }
    void operator+=(const Point& q) { x += q.x; y += q.y; }
    void operator-=(const Point& q) { x -= q.x; y -= q.y; }
    ll operator*(const Point& q) { return x * q.y - y * q.x; }
    ll prod_sc(const Point& q) const { return x * q.x + y * q.y; }
    bool operator==(const Point& q) { return x == q.x && y == q.y; }
    bool operator<(const Point& q) const { return x < q.x || (x == q.x && y < q.y); }
    ll dist2(const Point& q) const { return (x - q.x) * (x - q.x) + (y - q.y) * (y - q.y); }
    dd dist(const Point& q) const { return sqrt(dist2(q)); }
    // Check if a point lies on a line segment
    bool onSegment(Point p, Point q) const {
        ll aux=(p-*this)*(q-*this);
        if (aux==0&&min(p.x, q.x) <= x && x <= max(p.x, q.x) && min(p.y, q.y) <= y && y <= max(p.y, q.y)) {
            return true;
        }
        return false;
    }
    ll area(Point p1, Point p2) {
	return (p1.x -x) * (p2.y -y) -(p2.x -x) * (p1.y -y);
}

};

// Sign function to determine orientation
ll sgn(const ll& x) { 
    return x >= 0 ? (x ? 1LL : 0) : -1LL; 
}
// Check if two intervals overlap
bool inter(ll a, ll b, ll c, ll d) {
    if (a > b) swap(a, b);
    if (c > d) swap(c, d);
    return max(a, c) <= min(b, d);
}
// Check if two segments intersect
bool check_inter(Point& a, Point& b, Point& c, Point& d) {
    // Check if segments are collinear
    if ((a-c) * (d - c) == 0 && (b - c) * (d - c) == 0)
        return inter(a.x, b.x, c.x, d.x) && inter(a.y, b.y, c.y, d.y);
    // General case for segment intersection
    return (sgn((b - a) * (c - a)) != sgn((b - a) * (d - a)) &&
            sgn((d - c) * (a - c)) != sgn((d - c) * (b - c)));
}

vector<Point> monotone_chain(vector<Point>&points) {
    vector<Point> Hull;
	// sort with respect to the x and y coordinates
	sort(points.begin(), points.end());
	// distinct the points
	points.erase(unique(points.begin(), points.end()), points.end());
	ll n = points.size();

	// 1 or 2 points are always in the convex hull
	if (n < 3) {
		Hull = points;
		return Hull;
	}

	// lower hull
	for (ll i = 0; i < n; i++) {
		// if with the new point points[i], a right turn will be formed,
		// then we remove the last point in the hull and test further
		while (Hull.size() > 1 && 
        points[i].area(Hull[Hull.size() - 2],Hull.back()) <0) // we put < if we don't need the colinear point
		Hull.pop_back();                                       // else we put <= if we want them +colinear point
		// otherwise, add the point to the hull
		Hull.push_back(points[i]);
	}

	// upper hull, following the same logic as the lower hull
	auto lower_hull_length = Hull.size();
	for (ll i = n - 2; i >= 0; i--) {
		// we can only remove a point if there are still points left in the
		// upper hull
		while (Hull.size() > lower_hull_length && 
        points[i].area(Hull[Hull.size()-2],Hull.back()) < 0) // we put < if we don't need the colinear point
			Hull.pop_back();                                 // else we put <= if we want them +colinear point
		Hull.push_back(points[i]);
	}
	// delete point[0] that has been added twice
	Hull.pop_back();
    return Hull;
}


struct Point {
    dd x, y;
    Point(dd xx = 0.0, dd yy = 0.0) : x(xx), y(yy) {}
 
    bool operator<(const Point& other) const {
        if (x == other.x) return y < other.y;
        return x < other.x;
    }
};
 
dd distance(Point a, Point b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}
 
dd pointToSegmentDistance(Point p, Point a, Point b) {
    dd len = distance(a, b);
    if (len == 0.0) return distance(p, a);
    dd t = ((p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y)) / (len * len);
    if (t < 0.0) return distance(p, a);
    else if (t > 1.0) return distance(p, b);
    Point projection = {a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
    return distance(p, projection);
}
 
void affPO(Point P) {
    cout << P.x << " " << P.y << endl;
}
 

 pair<pair<dd, dd>, dd> findCircle(ll x1, ll y1, ll x2, ll y2, ll x3, ll y3) {
    ll x12 = x1 - x2;
    ll x13 = x1 - x3;
    ll y12 = y1 - y2;
    ll y13 = y1 - y3;
    ll y31 = y3 - y1;
    ll y21 = y2 - y1;
    ll x31 = x3 - x1;
    ll x21 = x2 - x1;

    ll sx13 = x1 * x1 - x3 * x3;
    ll sy13 = y1 * y1 - y3 * y3;
    ll sx21 = x2 * x2 - x1 * x1;
    ll sy21 = y2 * y2 - y1 * y1;

    dd f = ((sx13) * (x12) + (sy13) * (x12) + (sx21) * (x13) + (sy21) * (x13)) / (2.0 * ((y31) * (x12) - (y21) * (x13)));
    dd g = ((sx13) * (y12) + (sy13) * (y12) + (sx21) * (y13) + (sy21) * (y13)) / (2.0 * ((x31) * (y12) - (x21) * (y13)));

    dd c = -1.0 * x1 * x1 - y1 * y1 - 2.0 * g * x1 - 2.0 * f * y1;

    dd h = -g;
    dd k = -f;
    dd sqr_of_r = h * h + k * k - c;

    dd r = sqrt(sqr_of_r);
    return {{h, k}, r}; //// centre(h,k) rayon=r
    // dd pen1 = INF;
    // dd pen2 = INF;
    // if (x1 != x2) {
    //     pen1 = (dd)(y2 - y1) / (x2 - x1);
    // }                                          //// we put this in function sol
    // if (x2 != x3) {                             ///// EPS=1e-9
    //     pen2 = (dd)(y3 - y2) / (x3 - x2);
    // }
    // if (pen1 != pen2) {
    //     auto circle = findCircle(x1, y1, x2, y2, x3, y3);
    //     dd h = circle.ff.ff;
    //     dd k = circle.ff.ss;
    //     dd r = circle.ss;
    //     ll aux = 0;
    //     for (auto pr : vp) {
    //         dd tot = (h - pr.ff) * (h - pr.ff) + (k - pr.ss) * (k - pr.ss);
    //         if (abs(sqrt(tot) - r) < EPS) aux++;
    //     }
    //     ans = max(ans, aux);
    // }
}

// Circle-Line Intersection
vector<pair<dd, dd>> circle_line_intersection(dd cx, dd cy, dd r, dd a, dd b, dd c) {
    vector<pair<dd, dd>> intersections;
    dd x0 = -a * c / (a * a + b * b), y0 = -b * c / (a * a + b * b);
    if (c * c > r * r * (a * a + b * b) + EPS) {
        // No intersection
        return intersections;
    } else if (abs(c * c - r * r * (a * a + b * b)) < EPS) {
        // Tangent line
        intersections.push_back({x0 + cx, y0 + cy});
    } else {
        dd d = r * r - c * c / (a * a + b * b);
        dd mult = sqrt(d / (a * a + b * b));
        dd ax = x0 + b * mult, bx = x0 - b * mult;
        dd ay = y0 - a * mult, by = y0 + a * mult;
        intersections.push_back({ax + cx, ay + cy});
        intersections.push_back({bx + cx, by + cy});
    }
    return intersections;
}

// Circle-Circle Intersection
vector<pair<dd, dd>> circle_circle_intersection(dd x1, dd y1, dd r1, dd x2, dd y2, dd r2) {
    vector<pair<dd, dd>> intersections;
    dd d = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    if (d > r1 + r2 || d < abs(r1 - r2)) {
        // No intersection
        return intersections;
    }
    dd a = (r1 * r1 - r2 * r2 + d * d) / (2 * d);
    dd h = sqrt(r1 * r1 - a * a);
    dd x0 = x1 + a * (x2 - x1) / d;
    dd y0 = y1 + a * (y2 - y1) / d;
    dd rx = -(y2 - y1) * (h / d);
    dd ry = (x2 - x1) * (h / d);
    intersections.push_back({x0 + rx, y0 + ry});
    intersections.push_back({x0 - rx, y0 - ry});
    return intersections;
}



dd intersection_cercle_Area(dd X1,dd Y1,dd R1,dd X2,dd Y2,dd R2)
{
    dd Pi =acos(-1);
    dd d, alpha, beta, a1, a2;
    dd ans;
 
    // Calculate the euclidean distance
    // between the two points
    d = sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1));
 
    if (d > R1 + R2)
        ans = 0;
 
    else if (d <= (R1 - R2) && R1 >= R2)
        ans =Pi * R2 * R2;
 
    else if (d <= (R2 - R1) && R2 >= R1)
        ans =Pi * R1 * R1;
 
    else {
        alpha = acos((R1 * R1 + d * d - R2 * R2)
                     / (2 * R1 * d))
                * 2;
        beta = acos((R2 * R2 + d * d - R1 * R1)
                    / (2 * R2 * d))
               * 2;
        a1 = 0.5 * beta * R2 * R2
             - 0.5 * R2 * R2 * sin(beta);
        a2 = 0.5 * alpha * R1 * R1
             - 0.5 * R1 * R1 * sin(alpha);
        ans =a1 + a2;
    }
 
    return ans;
}


