#include <bits/stdc++.h>
/*
//#pragma GCC optimize("Ofast")
//#pragma GCC target("sse,sse2,sse3,ssse3,sse4")
*/
#define mp make_pair
#define mt make_tuple
#define all(x) (x).begin(), (x).end()
#define rall(x) (x).rbegin(), (x).rend()
#define sz(x) int((x).size())

using namespace std;

mt19937 gen_rand;

/*OUTPUT*/
/*PAIR*/
template<typename T, typename U>
ostream &operator<<(ostream &os, pair<T, U> &a) {
    os << "(";
    os << a.first << ", ";
    os << a.second;
    os << ")";
    return os;
}
/*PAIR*/

/*ARR*/
template<typename T>
ostream &operator<<(ostream &os, vector<T> &a) {
    os << "{";
    bool was = false;
    for (auto &x: a) {
        if (was) {
            os << ", ";
        }
        was = true;
        os << x;
    }
    os << "}";
    return os;
}
/*ARR*/
/*OUTPUT*/

/*DEBUG*/
template<typename T>
inline void debug(const char* sdbg, T x) {
    cerr << "The value of " << sdbg << " is " << x << "\n";
};

template<typename T, typename... U>
inline void debug(const char* sdbg, T h, U... t) {
    cerr << "The value of ";
    while (*sdbg != ',') {
        cerr << *sdbg++;
    }
    cerr << " is " << h << "\n";
    debug(sdbg + 1, t...);
    cerr << "\n";
}

#define value(...) debug(#__VA_ARGS__, __VA_ARGS__)
/*DEBUG*/

template<typename T>
T abs(T x) {
    if (x < 0) {
        return -x;
    } else {
        return x;
    }
}

template<typename T>
T sqr(T x) {
    return x * x;
};

const long long INF_FOR_SHORT_TIME = (long long)(1e9);

template<typename T>
T mul(T a, T b, T m) {
    if (a <= INF_FOR_SHORT_TIME && b <= INF_FOR_SHORT_TIME) {
        return (a * b) % m;
    } 
    T q = (long long)((long double)a * (long double)b / (long double)m);
    T r = a * b - q * m;
    return (r + 5 * m) % m;
}

template<typename T, typename U>
T binpow(T a, U n, T MOD) {
    T res = 1;
    while (n) {
        if (n & 1) {
            res = mul(res, a, MOD);
        } 
        a = mul(a, a, MOD);
        n >>= 1;
    }
    return res;
};

template<typename T>
int sign(T x) {
    if (x < 0) {
        return -1;
    } else if (x > 0) {
        return 1;
    } else {
        return 0;
    }
};

template<typename T>
T gcd(T a, T b) {
    if (a > b) {
        swap(a, b);
    }
    while (a) {
        b %= a;
        swap(a, b);
    }
    return b;
};

template<typename T>
bool uin(T &a, T b) {
    if (a > b) {
        a = b;
        return true;
    }
    return false;
};

template<typename T>
bool uax(T &a, T b) {
    if (a < b) {
        a = b;
        return true;
    }
    return false;
};

typedef unsigned int uint;
typedef unsigned long long ull;
typedef long long ll;
typedef long double ld;

/*curlink v1.9*/

int main() {
  return 0;
}
