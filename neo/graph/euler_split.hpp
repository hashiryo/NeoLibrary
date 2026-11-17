#pragma once
#include <array>
#include <cassert>
#include <vector>
// Euler 分割。euler_split(n, edges) は、頂点 n 個と無向辺 {u, v} の列を受け、辺ごとに 0 か 1 を返す。どの頂点でも、接する辺
// (自己ループは 2 本と数える) のうち 0 の本数と 1 の本数の差は、次数が奇数なら 1、偶数なら 0 になる。ただし、次数が全部偶数で
// 辺の数が奇数の連結成分では、その成分で番号がいちばん小さい頂点だけ差が 2 になる (三角形のように、どう塗っても避けられない)。
// 二部グラフにはそういう成分が無いので、どの頂点でも差は 1 以下になる。多重辺も受ける。
//
// 中身は、次数が奇数の頂点をすべて 1 つの補助の頂点とつないで次数を全部偶数にし、連結成分ごとに Euler 閉路をたどって、
// 通った順に辺を 0 と 1 に交互に塗る。閉路が頂点を通るたびに、入る辺と出る辺が違う色になるので、どの頂点でも 2 色の本数が
// 揃う。閉路の長さが奇数だと、始点でだけ同じ色が 2 本続く。補助の頂点のある成分は補助の頂点から始めるので、元の頂点は
// 揃ったままで、補助の辺を消すと次数が奇数の頂点だけ差が 1 になる。Euler 閉路は Hierholzer の方法を反復で書く。O(n + m)。
inline std::vector<int> euler_split(int n, const std::vector<std::array<int, 2>>& edges) {
 const int m= edges.size();
 std::vector<int> deg(n + 1, 0);
 for(auto [u, v]: edges) {
  assert(0 <= u && u < n && 0 <= v && v < n);
  ++deg[u], ++deg[v];
 }
 // 補助の辺 m + j は、次数が奇数の頂点 aux[j] と補助の頂点 n を結ぶ。
 std::vector<int> aux;
 for(int v= 0; v < n; ++v)
  if(deg[v] & 1) aux.push_back(v), ++deg[v], ++deg[n];
 const int M= m + aux.size();
 std::vector<int> off(n + 2, 0), adj(2 * (size_t)M);
 for(int v= 0; v <= n; ++v) off[v + 1]= off[v] + deg[v];
 std::vector<int> it(off.begin(), off.end() - 1);
 for(int e= 0; e < m; ++e) adj[it[edges[e][0]]++]= e, adj[it[edges[e][1]]++]= e;
 for(int j= 0; j < (int)aux.size(); ++j) adj[it[aux[j]]++]= m + j, adj[it[n]++]= m + j;
 it.assign(off.begin(), off.end() - 1);
 auto other= [&](int e, int x) { return e < m ? edges[e][0] ^ edges[e][1] ^ x : aux[e - m] ^ n ^ x; };
 std::vector<char> used(M, 0);
 std::vector<int> col(m, 0), stk;
 stk.reserve(2 * (size_t)M + 2);
 // s から Euler 閉路をたどり、辺を取り出した順に 0 と 1 を交互に塗る。stk は (頂点, 来た辺) を積む。取り出した辺の列は閉路を
 // 逆にたどったもので、続けて取り出した 2 本はその間の頂点で接し、最後と最初の 2 本は s で接する。
 auto run= [&](int s) {
  int k= 0;
  stk.push_back(s), stk.push_back(-1);
  while(!stk.empty()) {
   const int u= stk[stk.size() - 2];
   while(it[u] < off[u + 1] && used[adj[it[u]]]) ++it[u];
   if(it[u] == off[u + 1]) {
    const int e= stk.back();
    stk.pop_back(), stk.pop_back();
    if(e >= 0) {
     if(e < m) col[e]= k & 1;
     ++k;
    }
   } else {
    const int e= adj[it[u]++];
    used[e]= 1;
    stk.push_back(other(e, u)), stk.push_back(e);
   }
  }
 };
 if(!aux.empty()) run(n);
 for(int v= 0; v < n; ++v) {
  while(it[v] < off[v + 1] && used[adj[it[v]]]) ++it[v];
  if(it[v] < off[v + 1]) run(v);
 }
 return col;
}
