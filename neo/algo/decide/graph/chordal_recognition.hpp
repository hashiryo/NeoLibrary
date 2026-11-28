#pragma once
#include <algorithm>
#include <array>
#include <vector>
// 弦グラフ (長さ 4 以上の誘導閉路を持たない無向グラフ) の判定。弦グラフなら完全除去順序 (PEO) を、
// そうでなければ長さ 4 以上の誘導閉路を返す。最大基数探索 (MCS) を遅延バケットで O(N+M) に行い、
// PEO の検査は「x より後ろの隣接のうち最も早いもの」1 本との比較で済ませ、非弦のときは
// 証人 (x, z, y) から x と N(x)\{z} を避けた BFS で誘導閉路を作る。
// 入力は単純無向グラフ (自己ループと多重辺は未定義)。記録は algo-notes の notes/heavy-algorithms.md。
struct ChordalRecognition {
 bool chordal = true;
 std::vector<int> order;  // 弦グラフのときの PEO (長さ n)
 std::vector<int> cycle;  // 弦グラフでないときの誘導閉路 (巡回順、始点は繰り返さない)
};
inline ChordalRecognition chordal_recognition(int n, const std::vector<std::array<int, 2>>& edges) {
 ChordalRecognition res;
 std::vector<int> st(n + 1, 0), to(2 * edges.size());
 for(auto& [u, v]: edges) ++st[u + 1], ++st[v + 1];
 for(int i= 0; i < n; ++i) st[i + 1]+= st[i];
 {
  std::vector<int> p(st.begin(), st.end() - 1);
  for(auto& [u, v]: edges) to[p[u]++]= v, to[p[v]++]= u;
 }
 for(int i= 0; i < n; ++i) std::sort(to.begin() + st[i], to.begin() + st[i + 1]);
 auto adjacent= [&](int u, int v) { return std::binary_search(to.begin() + st[u], to.begin() + st[u + 1], v); };

 // MCS。重みごとのバケットに積み、取り出しで選択済みの頂点を捨てる。
 std::vector<int> weight(n, 0), pos(n, -1), ord;
 std::vector<std::vector<int>> buc(n);
 buc[0].reserve(n);
 for(int i= 0; i < n; ++i) buc[0].push_back(i);
 ord.reserve(n);
 int maxw= 0;
 for(int i= 0; i < n; ++i) {
  int x= -1;
  while(true) {
   while(maxw >= 0 && buc[maxw].empty()) --maxw;
   if(maxw < 0) break;
   x= buc[maxw].back(), buc[maxw].pop_back();
   if(pos[x] < 0) break;
  }
  pos[x]= i, ord.push_back(x);
  for(int k= st[x]; k < st[x + 1]; ++k)
   if(const int y= to[k]; pos[y] < 0) {
    if(++weight[y] > maxw) maxw= weight[y];
    buc[weight[y]].push_back(y);
   }
 }

 // x より前に選ばれた隣接のうち最も遅いもの z を取り、z より前の隣接 y が z に隣接して
 // いなければ (y, x, z) が証人になる。
 for(int x: ord) {
  int prev_pos= -1;
  for(int k= st[x]; k < st[x + 1]; ++k)
   if(const int y= to[k]; pos[y] < pos[x] && pos[y] > prev_pos) prev_pos= pos[y];
  if(prev_pos < 0) continue;
  const int z= ord[prev_pos];
  for(int k= st[x]; k < st[x + 1]; ++k) {
   const int y= to[k];
   if(pos[y] >= pos[z] || adjacent(z, y)) continue;
   // y から z へ、x と N(x)\{z} を避けた幅優先。見つけた道 z..y に y-x を足す。
   std::vector<char> banned(n, 0);
   banned[x]= 1;
   for(int j= st[x]; j < st[x + 1]; ++j)
    if(const int v= to[j]; v != z) banned[v]= 1;
   std::vector<int> from(n, -1), que;
   que.reserve(n);
   from[y]= y, que.push_back(y);
   for(std::size_t head= 0; head < que.size(); ++head) {
    const int u= que[head];
    if(u == z) break;
    for(int j= st[u]; j < st[u + 1]; ++j)
     if(const int v= to[j]; !banned[v] && from[v] < 0) from[v]= u, que.push_back(v);
   }
   if(from[z] < 0) continue;
   res.chordal= false;
   for(int u= z; u != y; u= from[u]) res.cycle.push_back(u);
   res.cycle.push_back(y), res.cycle.push_back(x);
   return res;
  }
 }
 res.order.assign(ord.rbegin(), ord.rend());
 return res;
}
