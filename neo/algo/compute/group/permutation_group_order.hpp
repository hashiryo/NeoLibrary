#pragma once
#include <cstdint>
#include <utility>
#include <vector>
// 置換 p (p[i] は i の行き先) が生成する群 G <= S_n の位数を mod で割った余りを返す。
// 合成は (p q)[i] = p[q[i]]。0 <= p[i] < n で、各 p は全単射。生成元に恒等置換や重複があってよい。
// incremental Schreier–Sims (Knuth の Algorithm A/B) で、各レベルの軌道の大きさの積を掛ける
// (割り算はしないので mod は素数でなくてよい)。最初に、2-サイクルがちょうど 1 つで他のサイクルが
// すべて奇数長の生成元から転置を作り、その共役でできる辺が全点を連結なら S_n として n! を返す。
// 記録は algo-notes の notes/heavy-algorithms.md。
inline unsigned permutation_group_order(int n, const std::vector<std::vector<int>>& generators, unsigned mod = 998244353) {
 if(mod <= 1) return 0;
 // 入力だけで S_n を証明する。失敗したら通常の鎖に落とす。
 auto certifies_symmetric= [&]() {
  if(n < 3 || generators.size() < 2) return false;
  int a= -1, b= -1;
  for(auto& g: generators) {
   std::vector<char> seen(n, 0);
   int two= 0, A= -1, B= -1;
   bool ok= true;
   for(int i= 0; i < n && ok; ++i) {
    if(seen[i]) continue;
    int len= 0, j= i;
    do {
     seen[j]= 1, j= g[j], ++len;
    } while(j != i);
    if(len == 2) ++two, A= i, B= g[i];
    else if(len % 2 == 0) ok= false;
   }
   if(ok && two == 1) {
    a= A, b= B;
    break;
   }
  }
  if(a < 0) return false;
  std::vector<int> par(n, -1);
  auto find= [&](int x) {
   int r= x;
   while(par[r] >= 0) r= par[r];
   while(x != r) {
    int p= par[x];
    par[x]= r, x= p;
   }
   return r;
  };
  auto unite= [&](int x, int y) {
   x= find(x), y= find(y);
   if(x == y) return false;
   if(par[x] > par[y]) std::swap(x, y);
   par[x]+= par[y], par[y]= x;
   return true;
  };
  unite(a, b);
  std::vector<std::pair<int, int>> edges{{a, b}};
  for(std::size_t ei= 0; ei < edges.size(); ++ei) {
   auto [u, v]= edges[ei];
   for(auto& p: generators) {
    const int x= p[u], y= p[v];
    if(!unite(x, y)) continue;
    edges.emplace_back(x, y);
    if((int)edges.size() == n - 1) return true;
   }
  }
  return false;
 };
 unsigned long long ans= 1;
 if(certifies_symmetric()) {
  for(int i= 2; i <= n; ++i) ans= ans * (unsigned)i % mod;
  return (unsigned)ans;
 }
 using Perm= std::vector<int>;
 // レベル i の代表元 T[i][x] (i を x に写す元)、その逆、軌道、生成元、処理済みの組の数。
 struct Chain {
  int n;
  std::vector<std::vector<Perm>> T, invT;
  std::vector<std::vector<char>> has;
  std::vector<std::vector<int>> orbit, done;
  std::vector<std::vector<Perm>> S;
  Chain(int n): n(n), T(n, std::vector<Perm>(n)), invT(n, std::vector<Perm>(n)), has(n, std::vector<char>(n, 0)), orbit(n), done(n, std::vector<int>(n, 0)), S(n) {
   Perm id(n);
   for(int i= 0; i < n; ++i) id[i]= i;
   for(int i= 0; i < n; ++i) T[i][i]= id, invT[i][i]= id, has[i][i]= 1, orbit[i].push_back(i);
  }
  static Perm compose(const Perm& a, const Perm& b) {
   const int n= a.size();
   Perm c(n);
   for(int i= 0; i < n; ++i) c[i]= a[b[i]];
   return c;
  }
  static Perm inverse(const Perm& p) {
   const int n= p.size();
   Perm q(n);
   for(int i= 0; i < n; ++i) q[p[i]]= i;
   return q;
  }
  bool contains(Perm g, int start) const {
   for(int i= start; i < n; ++i) {
    const int x= g[i];
    if(x == i) continue;
    if(!has[i][x]) return false;
    g= compose(invT[i][x], g);
   }
   return true;
  }
  void process(int i, Perm g) {
   if(i == n) return;
   const int x= g[i];
   if(!has[i][x]) {
    T[i][x]= g, invT[i][x]= inverse(g), has[i][x]= 1, orbit[i].push_back(x);
    return;
   }
   Perm h= (x == i) ? std::move(g) : compose(invT[i][x], g);
   if(!contains(h, i + 1)) extend(i + 1, std::move(h));
  }
  void extend(int i, Perm g) {
   S[i].push_back(std::move(g));
   for(std::size_t k= 0; k < orbit[i].size(); ++k) {
    const int x= orbit[i][k];
    while(done[i][x] < (int)S[i].size()) {
     const int j= done[i][x]++;
     process(i, compose(S[i][j], T[i][x]));
    }
   }
  }
 };
 Chain ch(n);
 for(auto& g: generators)
  if(!ch.contains(g, 0)) ch.extend(0, g);
 ans= 1;
 for(int i= 0; i < n; ++i) ans= ans * (unsigned)ch.orbit[i].size() % mod;
 return (unsigned)ans;
}
