#pragma once
#include <array>
#include <vector>
// 無向グラフを 2 色に塗り分ける。color[v] は 0 か 1 で、辺の両端は違う色になる。二部グラフでなければ (自己ループを
// 含めて) 空を返す。連結成分ごとに、番号の最も小さい頂点を 0 にする。BipartiteMatching の 2 つ目の入口に渡す塗り分けを
// 求めるのにも使う。
inline std::vector<int> bipartite_coloring(int n, const std::vector<std::array<int, 2>>& edges) {
 std::vector<int> st(n + 1), to(2 * edges.size()), col(n, -1), q(n);
 for(auto& [u, v]: edges) {
  if(u == v) return {};
  ++st[u + 1], ++st[v + 1];
 }
 for(int i= 0; i < n; ++i) st[i + 1]+= st[i];
 {
  std::vector<int> p(st.begin(), st.end() - 1);
  for(auto& [u, v]: edges) to[p[u]++]= v, to[p[v]++]= u;
 }
 for(int s= 0; s < n; ++s)
  if(col[s] < 0) {
   int h= 0, t= 0;
   col[s]= 0, q[t++]= s;
   while(h < t) {
    const int u= q[h++];
    for(int k= st[u]; k < st[u + 1]; ++k)
     if(const int w= to[k]; col[w] < 0) col[w]= col[u] ^ 1, q[t++]= w;
     else if(col[w] == col[u]) return {};
   }
  }
 return col;
}
