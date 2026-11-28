#pragma once
#include <algorithm>
#include <cstddef>
#include <vector>
// 順列 p (0 <= p[i] < n、重複なし) の common interval decomposition tree を作る。
// connected interval は [l, r) で値の集合も区間になるもの、strong interval はどの connected interval とも
// 交差しない (包含か非交) もので、strong interval 全体の包含のハッセ図がこの木になる。
// Uno–Yagiura の common interval 列挙を Xuan–Habib–Paul が木の構築に直したスタック法で、
// 「値域が左端より前を含むか」を値域の双方向リストで O(1) に見る。
// 節は構築順 (葉が 0..n-1、内部がそのあと) で、親は子より後の番号にもなる。問題文の p_i < i が
// 要るときは renumber_preorder を通す。記録は algo-notes の notes/heavy-algorithms.md。
struct PermutationTreeNode {
 int parent;  // 根は -1
 int l, r;    // 半開 [l, r)
 bool linear;  // linear 節なら true、prime 節なら false
};
inline std::vector<PermutationTreeNode> permutation_tree(const std::vector<int>& p) {
 const int n= p.size();
 std::vector<PermutationTreeNode> nodes;
 nodes.reserve(2 * n);
 std::vector<signed char> type;  // 0: prime、1: dec、2: inc、3: 葉
 type.reserve(2 * n);
 for(int i= 0; i < n; ++i) nodes.push_back({-1, i, i + 1, true}), type.push_back(3);
 if(n == 0) return nodes;
 // 値域の双方向リスト。後ろから抜くと、各値について位置 i より前での前後の値が残る。
 std::vector<int> lower(n), upper(n);
 for(int i= 0; i < n; ++i) lower[i]= i - 1, upper[i]= i + 1;
 for(int i= n - 1; i >= 0; --i) {
  const int v= p[i], lo= lower[v], hi= upper[v];
  if(lo >= 0) upper[lo]= hi;
  if(hi < n) lower[hi]= lo;
 }
 struct LeftBase {
  int l, vl, vr;
 };
 struct Common {
  int l, r;
 };
 std::vector<LeftBase> st;
 std::vector<Common> coms;
 st.reserve(n), coms.reserve(2 * n);
 for(int r= 1; r <= n; ++r) {
  int a= p[r - 1];
  LeftBase bs{r - 1, a, a + 1};
  while(!st.empty()) {
   if(bs.vl < st.back().vl) st.back().vl= bs.vl;
   if(bs.vr > st.back().vr) st.back().vr= bs.vr;
   LeftBase nx= st.back();
   if(!(lower[p[nx.l]] < nx.vl && nx.vr <= upper[p[nx.l]])) {
    // 値域 [vl, vr) が左端 l より前の値を含む。1 つ前の候補へ吸収する (一番下は l = 0 で来ない)。
    st.pop_back();
    LeftBase& nx2= st.back();
    if(nx.vl < nx2.vl) nx2.vl= nx.vl;
    if(nx.vr > nx2.vr) nx2.vr= nx.vr;
   } else if(nx.vr - nx.vl == r - nx.l) {
    bs= nx;
    st.pop_back();
    coms.push_back({nx.l, r});
   } else {
    break;
   }
  }
  st.push_back(bs);
 }
 while(st.size() >= 2) {
  LeftBase nx= st.back();
  st.pop_back();
  LeftBase& nx2= st.back();
  if(nx.vl < nx2.vl) nx2.vl= nx.vl;
  if(nx.vr > nx2.vr) nx2.vr= nx.vr;
  if(nx2.vr - nx2.vl == n - nx2.l) coms.push_back({nx2.l, n});
 }
 if(st.size() != 1) coms.push_back({0, n});
 // 確定した connected interval を左端の successor リストで組み立てる。
 std::vector<int> nodeid(n), sll(n);
 for(int i= 0; i < n; ++i) nodeid[i]= i, sll[i]= i + 1;
 for(const Common& com: coms) {
  const int m= sll[com.l];
  if(sll[m] == com.r) {
   // 2 タイル。左と同じ向きの linear なら右端を伸ばす。
   const int a= nodeid[com.l], b= nodeid[m];
   sll[com.l]= com.r;
   const signed char t= p[com.l] < p[com.r - 1] ? 2 : 1;
   if(type[a] == t) {
    nodes[b].parent= a, nodes[a].r= com.r;
   } else {
    const int c= nodes.size();
    nodes.push_back({-1, com.l, com.r, true}), type.push_back(t);
    nodes[a].parent= c, nodes[b].parent= c, nodeid[com.l]= c;
   }
  } else {
   // 3 タイル以上。全部を子にする prime。
   const int c= nodes.size();
   nodes.push_back({-1, com.l, com.r, false}), type.push_back(0);
   for(int x= com.l; x < com.r; x= sll[x]) nodes[nodeid[x]].parent= c;
   nodeid[com.l]= c, sll[com.l]= com.r;
  }
 }
 return nodes;
}
// 親が子より前になるように振り直す (行きがけ DFS、子は左端の昇順)。戻り値の節 i は parent < i を満たす。
inline std::vector<PermutationTreeNode> renumber_preorder(std::vector<PermutationTreeNode> nodes) {
 const int m= nodes.size();
 std::vector<std::vector<int>> ch(m);
 int root= -1;
 for(int i= 0; i < m; ++i)
  if(nodes[i].parent < 0) root= i;
  else ch[nodes[i].parent].push_back(i);
 for(auto& c: ch) std::sort(c.begin(), c.end(), [&](int a, int b) { return nodes[a].l < nodes[b].l; });
 std::vector<int> st{root}, order;
 order.reserve(m);
 while(!st.empty()) {
  const int x= st.back();
  st.pop_back(), order.push_back(x);
  for(auto it= ch[x].rbegin(); it != ch[x].rend(); ++it) st.push_back(*it);
 }
 std::vector<PermutationTreeNode> res(m);
 std::vector<int> newid(m);
 for(int i= 0; i < m; ++i) newid[order[i]]= i;
 for(int i= 0; i < m; ++i) {
  const PermutationTreeNode& v= nodes[order[i]];
  res[i]= {v.parent < 0 ? -1 : newid[v.parent], v.l, v.r, v.linear};
 }
 return res;
}
