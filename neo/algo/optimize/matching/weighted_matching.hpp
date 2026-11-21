#pragma once
#ifdef __x86_64__
#include <immintrin.h>
#else
#include <simde/x86/avx2.h>
#endif
#include <algorithm>
#include <array>
#include <cassert>
#include <limits>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
// 一般グラフの重み付きマッチング。辺は {u, v, w} の列で渡し、W は整数の型。多重辺と自己ループも受ける (自己ループは
// 使わない)。返すのは重みの和と、頂点ごとの相手 (いなければ -1)。重みの和が W に収まること。
// - max_weight_matching(n, edges): 重みの和が最大のマッチング (完全でなくてよい)。重みが 0 以下の辺は使わない。重みの最大を
//   W として 8W < 2^62 であること。
// - min_weight_perfect_matching(n, edges): 重みの和が最小の完全マッチング。無ければ nullopt。重みは負でもよい。重みの最大と
//   最小の差を D として 8 ((n / 2) D + 1) < 2^62 であること。
// 重み最大の完全マッチングと重み最小のマッチングは、重みの符号を反転して渡す。
// 頂点の数の 2 乗の行列を持つので、頂点は数千までにする。
//
// 中身は、よく知られた O(n^3) の主双対の花の方法 (頂点 2n の行列で花を持つ形) に、交互森を増加路ごとに作り直さない手を
// 足したもの。増加路を流したら、その 2 本の木の印だけを消し、ほかの木の印と slack はそのまま使い続ける。双対は動かした量の
// 和 T を使って遅延で持ち、slack も縮約費用に T の倍数を足した値で持つので、双対を動かす手間は T を足すだけになる。木を
// 壊すと、壊した木の頂点を指す slack が古くなるが、古い値は本当の値以下なので、双対の動かし方を決めるときに最小を与えた
// ものだけを確かめて付け直す。元の頂点についての走査 (偶の頂点の行、slack の付け直し、双対の動かし方の最小、動かした
// あとのきつい slack の検出、花の行の計算) は AVX2 で回す。重みの最大 W が 8W < 2^30 なら、値を int32 で持って 8 列ずつ、
// そうでなければ i64 で 4 列ずつ回す。
// 重み最小の完全マッチングは、重みの最大と最小を hi と lo、C = (n / 2) (hi - lo) + lo + 1 として、重みを C - w に替えた重み
// 最大のマッチングで解く。どの完全マッチングも、どの完全でないマッチングより重くなる。
// 方式は procon-judge の yosupo-general-weighted-matching で比べて選んだ (persist_i32 と persist_avx2)。記録は algo-notes の
// notes/network-flow.md。
template <class W> struct WeightedMatching {
 W weight;               // 重みの和
 std::vector<int> mate;  // 頂点ごとの相手。いなければ -1
};
namespace weighted_matching_internal {
// 頂点は内部で 1 から n。重みが正の辺だけを持つ (重み 0 は辺が無いことを表す)。add_edge は多重辺のうち重い方を残す。
// solve() は重みの和が最大のマッチングの、頂点 (0 始まり) ごとの相手を返す。Solver64 は値を i64、Solver32 は int32 で持つ。
struct Solver64 {
 using i64= long long;
 static constexpr i64 INF= (i64)1 << 62;
 struct E {
  int u, v;
 };
 int n, nx, N2, n1;
 std::vector<E> g;  // n × (2n + 1)。花 b の行 b - n - 1 に、花の中の端と相手の端を持つ
 std::vector<i64> wv, lb, key, bd;
 std::vector<int> sv, match, slack, st, pa, bp, ftime, S, vis, q, rt, cand, dead, osince, kst;
 std::vector<std::vector<int>> flo;
 int qh= 0, tvis= 0, clk= 0, fclk= 0;
 i64 T= 0;
 E& row(int b, int x) { return g[(size_t)(b - n - 1) * N2 + x]; }
 const E& row(int b, int x) const { return g[(size_t)(b - n - 1) * N2 + x]; }
 // a と b を結ぶ代表の辺 (u が a の側)。あとから作った花の行が新しい。
 E rep(int a, int b) const {
  if(a <= n && b <= n) return E{a, b};
  if(a > n && (b <= n || ftime[a] > ftime[b])) return row(a, b);
  const E& r= row(b, a);
  return E{r.v, r.u};
 }
 // 花 b の子のうち、元の頂点 m を含むもの。
 int child_of(int b, int m) const {
  while(bp[m] != b) m= bp[m];
  return m;
 }
 i64 w0(int u, int v) const { return wv[u * n1 + v]; }
 i64 lab(int x) const { return lb[x] + sv[x] * T; }
 i64 dlt(const E& x) const { return lab(x.u) + lab(x.v) - w0(x.u, x.v) * 2; }
 i64 cf(int x) const { return S[x] == 0 ? 2 : 1; }
 explicit Solver64(int n): n(n), nx(n), N2(2 * n + 1), n1(n + 1), g((size_t)n * N2), wv((size_t)n1 * n1), lb(N2), key(N2, INF), bd(N2), sv(N2), match(N2), slack(N2), st(N2), pa(N2), bp(N2), ftime(N2), S(N2, -1), vis(N2), rt(N2), osince(n1), kst(N2), flo(N2) {}
 void add_edge(int u, int v, i64 w) {  // 0 始まり
  ++u, ++v;
  if(w > w0(u, v)) wv[u * n1 + v]= wv[v * n1 + u]= w;
 }
 // 頂点の印を変えるときに、中の頂点の双対の持ち方を付け直す。c は頂点の係数 (-1, 0, +1)。偶になった時刻を osince に残す。
 void set_vertex_coef(int x, int c) {
  if(x <= n) {
   if(c == -1 && sv[x] != -1) osince[x]= ++clk;
   lb[x]= lab(x) - c * T, sv[x]= c;
   return;
  }
  for(int y: flo[x]) set_vertex_coef(y, c);
 }
 // 頂点の花 b の双対の係数を付け直す (c は -2, 0, +2)。
 void set_blossom_coef(int b, int c) {
  if(b > n) lb[b]= lab(b) - c * T, sv[b]= c;
 }
 // 頂点 (または花) x の印を s にする。双対の係数も付け直す。
 void set_label(int x, int s) {
  S[x]= s;
  set_vertex_coef(x, s == 0 ? -1 : s == 1 ? 1 : 0);
  set_blossom_coef(x, s == 0 ? 2 : s == 1 ? -2 : 0);
 }
 // slack[x] が、記録したときからずっと偶のままの頂点を指しているか。途中で偶でなくなった頂点は、縮約費用が想定どおりに
 // 減っていないので、偶に戻っていても古い。
 bool slack_valid(int x) const { return S[st[slack[x]]] == 0 && st[slack[x]] != x && osince[slack[x]] <= kst[x]; }
 // slack は元の頂点と花のどちらも、縮約費用に cf(x) T を足した値 key[x] で持つ。偶のままの頂点からの値は変わらない。
 void reset_slack(int x) { slack[x]= 0, key[x]= INF; }
 i64 edge_d(int u, int x) const { return x <= n ? lab(u) + lab(x) - 2 * w0(u, x) : dlt(row(x, u)); }
 // 古い slack と比べて小さければ書き換えてよい。古い slack の値は、それに負けた偶のままの頂点のどれの値よりも小さいので、
 // それより小さい値は本当の最小になる。大きければ捨てる (古い slack は、双対の動かし方を決めるときに付け直す)。
 void update_slack(int u, int x) {
  if(const i64 c= edge_d(u, x) + cf(x) * T; c < key[x]) key[x]= c, slack[x]= u, kst[x]= clk;
 }
 static __m256i ld(const i64* p) { return _mm256_loadu_si256((const __m256i*)p); }
 static __m256i ld32(const int* p) { return _mm256_cvtepi32_epi64(_mm_loadu_si128((const __m128i*)p)); }
 // 4 列の最小とその列を、スカラーの値と合わせて返す。
 static void reduce_min(__m256i bv, __m256i bi, i64& b, int& bx) {
  alignas(32) i64 v[4], id[4];
  _mm256_store_si256((__m256i*)v, bv), _mm256_store_si256((__m256i*)id, bi);
  for(int j= 0; j < 4; ++j)
   if(v[j] < b) b= v[j], bx= (int)id[j];
 }
 void set_slack(int x) {
  reset_slack(x);
  if(x <= n) {
   // 偶の頂点 u は sv[u] == -1 で、lab(u) = lb[u] - T。x と同じ頂点に入る u は x 自身だけで、自己ループは無い。
   const i64* wr= &wv[x * n1];
   const i64 base= lab(x) - T + cf(x) * T;
   const __m256i vb= _mm256_set1_epi64x(base), zero= _mm256_setzero_si256(), m1= _mm256_set1_epi64x(-1);
   __m256i bv= _mm256_set1_epi64x(INF), bi= zero, id= _mm256_setr_epi64x(1, 2, 3, 4);
   const __m256i four= _mm256_set1_epi64x(4);
   int u= 1;
   for(; u + 3 <= n; u+= 4, id= _mm256_add_epi64(id, four)) {
    const __m256i w= ld(wr + u);
    const __m256i act= _mm256_and_si256(_mm256_cmpgt_epi64(w, zero), _mm256_cmpeq_epi64(ld32(sv.data() + u), m1));
    const __m256i val= _mm256_sub_epi64(_mm256_add_epi64(ld(lb.data() + u), vb), _mm256_add_epi64(w, w));
    const __m256i better= _mm256_and_si256(act, _mm256_cmpgt_epi64(bv, val));
    bv= _mm256_blendv_epi8(bv, val, better), bi= _mm256_blendv_epi8(bi, id, better);
   }
   i64 b= INF;
   int bx= 0;
   reduce_min(bv, bi, b, bx);
   for(; u <= n; ++u)
    if(wr[u] > 0 && sv[u] == -1)
     if(const i64 val= lb[u] + base - 2 * wr[u]; val < b) b= val, bx= u;
   if(bx) key[x]= b, slack[x]= bx, kst[x]= clk;
  } else
   for(int u= 1; u <= n; ++u)
    if(const E& ex= row(x, u); w0(ex.u, ex.v) > 0 && st[u] != x && S[st[u]] == 0) update_slack(u, x);
 }
 i64 slack_value(int x) const { return key[x] - cf(x) * T; }
 void q_push(int x) {
  if(x <= n) q.push_back(x);
  else
   for(int y: flo[x]) q_push(y);
 }
 void set_st(int x, int b) {
  st[x]= b;
  if(x > n)
   for(int y: flo[x]) set_st(y, b);
 }
 int get_pr(int b, int xr) {
  int pr= std::find(flo[b].begin(), flo[b].end(), xr) - flo[b].begin();
  if(pr % 2 == 1) {
   std::reverse(flo[b].begin() + 1, flo[b].end());
   return (int)flo[b].size() - pr;
  }
  return pr;
 }
 void set_match(int u, int v) {
  const E r= rep(u, v);
  match[u]= r.v;
  if(u <= n) return;
  const int xr= child_of(u, r.u), pr= get_pr(u, xr);
  for(int i= 0; i < pr; ++i) set_match(flo[u][i], flo[u][i ^ 1]);
  set_match(xr, v);
  std::rotate(flo[u].begin(), flo[u].begin() + pr, flo[u].end());
 }
 void augment(int u, int v) {
  for(;;) {
   const int xnv= st[match[u]];
   set_match(u, v);
   if(!xnv) return;
   set_match(xnv, st[pa[xnv]]);
   u= st[pa[xnv]], v= xnv;
  }
 }
 int get_lca(int u, int v) {
  for(++tvis; u || v; std::swap(u, v)) {
   if(u == 0) continue;
   if(vis[u] == tvis) return u;
   vis[u]= tvis;
   u= st[match[u]];
   if(u) u= st[pa[u]];
  }
  return 0;
 }
 void add_blossom(int u, int lca, int v) {
  int b= n + 1;
  while(b <= nx && st[b]) ++b;
  if(b > nx) ++nx;
  match[b]= match[lca];
  flo[b].clear();
  flo[b].push_back(lca);
  for(int x= u, y; x != lca; x= st[pa[y]]) flo[b].push_back(x), flo[b].push_back(y= st[match[x]]), q_push(y);
  std::reverse(flo[b].begin() + 1, flo[b].end());
  for(int x= v, y; x != lca; x= st[pa[y]]) flo[b].push_back(x), flo[b].push_back(y= st[match[x]]), q_push(y);
  // 中に入る花の双対は止め、中の頂点はすべて偶にする。新しい花の双対は 0 から始める。
  for(int x: flo[b]) set_blossom_coef(x, 0), set_vertex_coef(x, -1), S[x]= -1, bp[x]= b;
  rt[b]= rt[lca], S[b]= 0, sv[b]= 2, lb[b]= -2 * T, bp[b]= 0;
  set_st(b, b);
  // 花の行。どの子の頂点も偶なので、元の頂点 x への辺 (m, x) の縮約費用は lb[m] - T + lab(x) - 2w(m, x)。
  E* rb= &row(b, 0);
  std::fill(rb, rb + N2, E{0, 0});
  std::fill(bd.begin(), bd.begin() + nx + 1, INF);
  const __m256i vT= _mm256_set1_epi64x(T), zero= _mm256_setzero_si256(), one64= _mm256_set1_epi64x(1), m1= _mm256_set1_epi64x(-1);
  const __m256i vn1= _mm256_set1_epi64x(n1), lo32= _mm256_set1_epi64x(0xffffffffLL);
  for(int xs: flo[b]) {
   const E* rx= xs > n ? &row(xs, 0) : nullptr;
   const __m256i vxs= _mm256_set1_epi64x(xs), lbxs= _mm256_set1_epi64x(xs <= n ? lb[xs] : 0);
   int x= 1;
   for(; x + 3 <= n; x+= 4) {
    const __m256i id= _mm256_add_epi64(_mm256_set1_epi64x(x), _mm256_setr_epi64x(0, 1, 2, 3));
    __m256i m, w, lm;
    if(!rx) m= vxs, w= ld(wv.data() + (size_t)xs * n1 + x), lm= lbxs;
    else {
     m= _mm256_and_si256(_mm256_loadu_si256((const __m256i*)(rx + x)), lo32);
     w= _mm256_i64gather_epi64((const long long*)wv.data(), _mm256_add_epi64(_mm256_mul_epu32(m, vn1), id), 8);
     lm= _mm256_i64gather_epi64((const long long*)lb.data(), m, 8);
    }
    const __m256i s4= ld32(sv.data() + x);
    const __m256i labx= _mm256_sub_epi64(_mm256_add_epi64(ld(lb.data() + x), _mm256_and_si256(_mm256_cmpeq_epi64(s4, one64), vT)), _mm256_and_si256(_mm256_cmpeq_epi64(s4, m1), vT));
    const __m256i d= _mm256_sub_epi64(_mm256_add_epi64(_mm256_sub_epi64(lm, vT), labx), _mm256_add_epi64(w, w));
    const __m256i cur= ld(bd.data() + x);
    const __m256i better= _mm256_and_si256(_mm256_cmpgt_epi64(w, zero), _mm256_cmpgt_epi64(cur, d));
    if(_mm256_testz_si256(better, better)) continue;
    _mm256_storeu_si256((__m256i*)(bd.data() + x), _mm256_blendv_epi8(cur, d, better));
    const __m256i ev= _mm256_or_si256(m, _mm256_slli_epi64(id, 32));
    const __m256i old= _mm256_loadu_si256((const __m256i*)(rb + x));
    _mm256_storeu_si256((__m256i*)(rb + x), _mm256_blendv_epi8(old, ev, better));
   }
   for(; x <= nx; ++x) {
    if(x == b) continue;
    const E c= rep(xs, x);
    if(!c.u || w0(c.u, c.v) == 0) continue;
    if(const i64 d= dlt(c); d < bd[x]) bd[x]= d, rb[x]= c;
   }
  }
  ftime[b]= ++fclk;
  set_slack(b);
 }
 void expand_blossom(int b) {
  for(int y: flo[b]) set_st(y, y);
  const int xr= child_of(b, rep(b, pa[b]).u), pr= get_pr(b, xr);
  for(int y: flo[b]) bp[y]= 0;
  for(int i= 0; i < pr; i+= 2) {
   const int xs= flo[b][i], xns= flo[b][i + 1];
   pa[xs]= rep(xns, xs).u;
   set_label(xs, 1), set_label(xns, 0), rt[xs]= rt[xns]= rt[b];
   reset_slack(xs), set_slack(xns);
   q_push(xns);
  }
  set_label(xr, 1), rt[xr]= rt[b], pa[xr]= pa[b];
  reset_slack(xr);
  for(size_t i= pr + 1; i < flo[b].size(); ++i) {
   const int xs= flo[b][i];
   set_label(xs, -1), set_slack(xs);
  }
  st[b]= 0, S[b]= -1;
 }
 // 増加路を流したら、両端の木の根を dead に入れて true を返す。
 bool on_found_edge(const E& x) {
  const int u= st[x.u], v= st[x.v];
  if(S[v] == -1) {
   pa[v]= x.u, set_label(v, 1), rt[v]= rt[u];
   const int nu= st[match[v]];
   reset_slack(v), reset_slack(nu);
   set_label(nu, 0), rt[nu]= rt[u], q_push(nu);
  } else if(S[v] == 0) {
   const int lca= get_lca(u, v);
   if(!lca) {
    dead.push_back(rt[u]), dead.push_back(rt[v]);
    return augment(u, v), augment(v, u), true;
   }
   add_blossom(u, lca, v);
  }
  return false;
 }
 // 流した 2 本の木を壊す。印を消してから、生きている木の偶の頂点からの slack を付け直し、きつい辺があれば伸ばす。
 void dissolve() {
  const int ra= dead[0], rb= dead[1];
  dead.clear();
  cand.clear();
  for(int x= 1; x <= nx; ++x)
   if(st[x] == x && S[x] != -1 && (rt[x] == ra || rt[x] == rb)) set_label(x, -1), reset_slack(x), cand.push_back(x);
  for(int x: cand) {
   if(S[x] != -1) continue;  // 先の頂点から伸びた木に拾われた
   set_slack(x);
   if(slack[x] && slack_value(x) == 0) on_found_edge(rep(slack[x], x));
  }
 }
 // 偶の頂点 u の行を走査する。増加路を流したら true。
 bool scan(int u) {
  const i64* wr= &wv[u * n1];
  const i64 lu= lab(u);
  auto one= [&](int v) {
   if(wr[v] > 0 && st[v] == v && st[u] != v && S[v] != 1) {
    if(lu + lab(v) - 2 * wr[v] == 0) return on_found_edge(E{u, v});
    update_slack(u, v);
   }
   return false;
  };
  const __m256i vlu= _mm256_set1_epi64x(lu), vT= _mm256_set1_epi64x(T), zero= _mm256_setzero_si256();
  const __m256i one64= _mm256_set1_epi64x(1), m1= _mm256_set1_epi64x(-1);
  int v= 1;
  for(; v + 3 <= n; v+= 4) {
   const __m256i w= ld(wr + v), s4= ld32(sv.data() + v), S4= ld32(S.data() + v);
   const __m256i id= _mm256_add_epi64(_mm256_set1_epi64x(v), _mm256_setr_epi64x(0, 1, 2, 3));
   const __m256i labv= _mm256_sub_epi64(_mm256_add_epi64(ld(lb.data() + v), _mm256_and_si256(_mm256_cmpeq_epi64(s4, one64), vT)), _mm256_and_si256(_mm256_cmpeq_epi64(s4, m1), vT));
   const __m256i dd= _mm256_sub_epi64(_mm256_add_epi64(vlu, labv), _mm256_add_epi64(w, w));
   const __m256i act= _mm256_andnot_si256(_mm256_cmpeq_epi64(S4, one64), _mm256_and_si256(_mm256_cmpgt_epi64(w, zero), _mm256_cmpeq_epi64(ld32(st.data() + v), id)));
   if(!_mm256_testz_si256(act, _mm256_cmpeq_epi64(dd, zero))) {
    for(int j= 0; j < 4; ++j)
     if(one(v + j)) return true;
    continue;
   }
   const __m256i cand= _mm256_add_epi64(_mm256_add_epi64(dd, vT), _mm256_and_si256(_mm256_cmpeq_epi64(S4, zero), vT));
   const __m256i kv= ld(key.data() + v);
   const __m256i upd= _mm256_and_si256(act, _mm256_cmpgt_epi64(kv, cand));
   if(_mm256_testz_si256(upd, upd)) continue;
   _mm256_storeu_si256((__m256i*)(key.data() + v), _mm256_blendv_epi8(kv, cand, upd));
   for(int m= _mm256_movemask_pd(_mm256_castsi256_pd(upd)); m; m&= m - 1) slack[v + __builtin_ctz(m)]= u, kst[v + __builtin_ctz(m)]= clk;
  }
  for(; v <= n; ++v)
   if(one(v)) return true;
  for(int b= n + 1; b <= nx; ++b)
   if(st[b] == b && S[b] != 1 && st[u] != b && row(b, u).u) {
    if(dlt(row(b, u)) == 0) {
     if(on_found_edge(rep(u, b))) return true;
    } else update_slack(u, b);
   }
  return false;
 }
 std::vector<int> solve() {
  nx= n;
  for(int u= 0; u <= n; ++u) st[u]= u, flo[u].clear();
  i64 wmax= 0;
  for(int u= 1; u <= n; ++u)
   for(int v= 1; v <= n; ++v) wmax= std::max(wmax, w0(u, v));
  for(int u= 1; u <= n; ++u) lb[u]= wmax, rt[u]= u, pa[u]= 0, set_label(u, 0), q.push_back(u);
  for(;;) {
   while(qh < (int)q.size()) {
    const int u= q[qh++];
    if(S[st[u]] != 0) continue;
    if(scan(u)) dissolve();
   }
   // 双対の動かし方。slack は 2 倍の尺度で比べる (印の無い頂点への辺は縮約費用の 2 倍、偶どうしは縮約費用、奇の花は双対)。
   i64 best= INF;
   // 元の頂点の最小。最小を与えたものが古ければ付け直してやり直す (古い slack は本当の値以下)。
   for(;;) {
    const __m256i vT= _mm256_set1_epi64x(T), zero= _mm256_setzero_si256(), one64= _mm256_set1_epi64x(1), m1= _mm256_set1_epi64x(-1);
    const __m256i vinf= _mm256_set1_epi64x(INF), four= _mm256_set1_epi64x(4);
    __m256i bv= vinf, bi= zero, id= _mm256_setr_epi64x(1, 2, 3, 4);
    int x= 1;
    for(; x + 3 <= n; x+= 4, id= _mm256_add_epi64(id, four)) {
     const __m256i S4= ld32(S.data() + x), kv= ld(key.data() + x);
     const __m256i act= _mm256_andnot_si256(_mm256_cmpeq_epi64(S4, one64), _mm256_and_si256(_mm256_cmpeq_epi64(ld32(st.data() + x), id), _mm256_cmpgt_epi64(vinf, kv)));
     const __m256i val= _mm256_sub_epi64(_mm256_sub_epi64(kv, vT), _mm256_and_si256(_mm256_cmpeq_epi64(S4, zero), vT));
     const __m256i sc= _mm256_add_epi64(val, _mm256_and_si256(_mm256_cmpeq_epi64(S4, m1), val));
     const __m256i better= _mm256_and_si256(act, _mm256_cmpgt_epi64(bv, sc));
     bv= _mm256_blendv_epi8(bv, sc, better), bi= _mm256_blendv_epi8(bi, id, better);
    }
    i64 b= INF;
    int bx= 0;
    reduce_min(bv, bi, b, bx);
    for(; x <= n; ++x)
     if(st[x] == x && S[x] != 1 && slack[x])
      if(const i64 val= slack_value(x) * (S[x] == -1 ? 2 : 1); val < b) b= val, bx= x;
    if(bx && !slack_valid(bx)) {
     set_slack(bx);
     continue;
    }
    best= b;
    break;
   }
   for(int x= n + 1; x <= nx; ++x) {
    if(st[x] != x || S[x] == 1 || !slack[x]) continue;
    i64 val= slack_value(x) * (S[x] == -1 ? 2 : 1);
    if(val >= best) continue;
    if(!slack_valid(x)) {
     set_slack(x);
     if(!slack[x]) continue;
     val= slack_value(x) * (S[x] == -1 ? 2 : 1);
     if(val >= best) continue;
    }
    best= val;
   }
   for(int b= n + 1; b <= nx; ++b)
    if(st[b] == b && S[b] == 1) best= std::min(best, lab(b));
   const i64 d= best / 2;
   // 偶の頂点の双対の最小 (lb - T)。これが d 以下なら終わり。偶の頂点が無ければ空いている頂点も無い。
   {
    const __m256i m1= _mm256_set1_epi64x(-1);
    __m256i mv= _mm256_set1_epi64x(INF);
    int u= 1;
    for(; u + 3 <= n; u+= 4) {
     const __m256i l= ld(lb.data() + u);
     const __m256i act= _mm256_and_si256(_mm256_cmpeq_epi64(ld32(sv.data() + u), m1), _mm256_cmpgt_epi64(mv, l));
     mv= _mm256_blendv_epi8(mv, l, act);
    }
    i64 mlb= INF;
    int dummy= 0;
    reduce_min(mv, _mm256_setzero_si256(), mlb, dummy);
    for(; u <= n; ++u)
     if(sv[u] == -1) mlb= std::min(mlb, lb[u]);
    if(mlb >= INF || best >= INF || mlb - T <= d) break;
   }
   T+= d;
   {
    const __m256i vT= _mm256_set1_epi64x(T), zero= _mm256_setzero_si256(), one64= _mm256_set1_epi64x(1);
    int x= 1;
    auto one= [&](int y) {
     if(st[y] == y && S[y] != 1 && slack[y] && slack_value(y) == 0 && slack_valid(y))
      if(on_found_edge(rep(slack[y], y))) dissolve();
    };
    for(; x + 3 <= n; x+= 4) {
     const __m256i S4= ld32(S.data() + x);
     const __m256i val= _mm256_sub_epi64(_mm256_sub_epi64(ld(key.data() + x), vT), _mm256_and_si256(_mm256_cmpeq_epi64(S4, zero), vT));
     const __m256i hit= _mm256_andnot_si256(_mm256_cmpeq_epi64(S4, one64), _mm256_cmpeq_epi64(val, zero));
     for(int m= _mm256_movemask_pd(_mm256_castsi256_pd(hit)); m; m&= m - 1) one(x + __builtin_ctz(m));
    }
    for(; x <= nx; ++x) one(x);
   }
   for(int b= n + 1; b <= nx; ++b)
    if(st[b] == b && S[b] == 1 && lab(b) == 0) expand_blossom(b);
  }
  std::vector<int> mt(n, -1);
  for(int u= 1; u <= n; ++u)
   if(match[u]) mt[u - 1]= match[u] - 1;
  return mt;
 }
};
struct Solver32 {
 using V= int;
 static constexpr V INF= 1 << 30;
 struct E {
  int u, v;
 };
 int n, nx, N2, n1;
 std::vector<E> g;  // n × (2n + 1)。花 b の行 b - n - 1 に、花の中の端と相手の端を持つ
 std::vector<V> wv, lb, key, bd;
 std::vector<int> sv, match, slack, st, pa, bp, ftime, S, vis, q, rt, cand, dead, osince, kst;
 std::vector<std::vector<int>> flo;
 int qh= 0, tvis= 0, clk= 0, fclk= 0;
 V T= 0;
 E& row(int b, int x) { return g[(size_t)(b - n - 1) * N2 + x]; }
 const E& row(int b, int x) const { return g[(size_t)(b - n - 1) * N2 + x]; }
 // a と b を結ぶ代表の辺 (u が a の側)。あとから作った花の行が新しい。
 E rep(int a, int b) const {
  if(a <= n && b <= n) return E{a, b};
  if(a > n && (b <= n || ftime[a] > ftime[b])) return row(a, b);
  const E& r= row(b, a);
  return E{r.v, r.u};
 }
 // 花 b の子のうち、元の頂点 m を含むもの。
 int child_of(int b, int m) const {
  while(bp[m] != b) m= bp[m];
  return m;
 }
 V w0(int u, int v) const { return wv[u * n1 + v]; }
 V lab(int x) const { return lb[x] + sv[x] * T; }
 V dlt(const E& x) const { return lab(x.u) + lab(x.v) - w0(x.u, x.v) * 2; }
 V cf(int x) const { return S[x] == 0 ? 2 : 1; }
 explicit Solver32(int n): n(n), nx(n), N2(2 * n + 1), n1(n + 1), g((size_t)n * N2), wv((size_t)n1 * n1), lb(N2), key(N2, INF), bd(N2), sv(N2), match(N2), slack(N2), st(N2), pa(N2), bp(N2), ftime(N2), S(N2, -1), vis(N2), rt(N2), osince(n1), kst(N2), flo(N2) {}
 void add_edge(int u, int v, V w) {  // 0 始まり
  ++u, ++v;
  if(w > w0(u, v)) wv[u * n1 + v]= wv[v * n1 + u]= w;
 }
 // 頂点の印を変えるときに、中の頂点の双対の持ち方を付け直す。c は頂点の係数 (-1, 0, +1)。偶になった時刻を osince に残す。
 void set_vertex_coef(int x, int c) {
  if(x <= n) {
   if(c == -1 && sv[x] != -1) osince[x]= ++clk;
   lb[x]= lab(x) - c * T, sv[x]= c;
   return;
  }
  for(int y: flo[x]) set_vertex_coef(y, c);
 }
 // 頂点の花 b の双対の係数を付け直す (c は -2, 0, +2)。
 void set_blossom_coef(int b, int c) {
  if(b > n) lb[b]= lab(b) - c * T, sv[b]= c;
 }
 // 頂点 (または花) x の印を s にする。双対の係数も付け直す。
 void set_label(int x, int s) {
  S[x]= s;
  set_vertex_coef(x, s == 0 ? -1 : s == 1 ? 1 : 0);
  set_blossom_coef(x, s == 0 ? 2 : s == 1 ? -2 : 0);
 }
 // slack[x] が、記録したときからずっと偶のままの頂点を指しているか。途中で偶でなくなった頂点は、縮約費用が想定どおりに
 // 減っていないので、偶に戻っていても古い。
 bool slack_valid(int x) const { return S[st[slack[x]]] == 0 && st[slack[x]] != x && osince[slack[x]] <= kst[x]; }
 // slack は元の頂点と花のどちらも、縮約費用に cf(x) T を足した値 key[x] で持つ。偶のままの頂点からの値は変わらない。
 void reset_slack(int x) { slack[x]= 0, key[x]= INF; }
 V edge_d(int u, int x) const { return x <= n ? lab(u) + lab(x) - 2 * w0(u, x) : dlt(row(x, u)); }
 // 古い slack と比べて小さければ書き換えてよい。古い slack の値は、それに負けた偶のままの頂点のどれの値よりも小さいので、
 // それより小さい値は本当の最小になる。大きければ捨てる (古い slack は、双対の動かし方を決めるときに付け直す)。
 void update_slack(int u, int x) {
  if(const V c= edge_d(u, x) + cf(x) * T; c < key[x]) key[x]= c, slack[x]= u, kst[x]= clk;
 }
 static __m256i ld(const int* p) { return _mm256_loadu_si256((const __m256i*)p); }
 // 4 列の最小とその列を、スカラーの値と合わせて返す。
 static void reduce_min(__m256i bv, __m256i bi, V& b, int& bx) {
  alignas(32) int v[8], id[8];
  _mm256_store_si256((__m256i*)v, bv), _mm256_store_si256((__m256i*)id, bi);
  for(int j= 0; j < 8; ++j)
   if(v[j] < b) b= v[j], bx= (int)id[j];
 }
 void set_slack(int x) {
  reset_slack(x);
  if(x <= n) {
   // 偶の頂点 u は sv[u] == -1 で、lab(u) = lb[u] - T。x と同じ頂点に入る u は x 自身だけで、自己ループは無い。
   const V* wr= &wv[x * n1];
   const V base= lab(x) - T + cf(x) * T;
   const __m256i vb= _mm256_set1_epi32(base), zero= _mm256_setzero_si256(), m1= _mm256_set1_epi32(-1);
   __m256i bv= _mm256_set1_epi32(INF), bi= zero, id= _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8);
   const __m256i eight= _mm256_set1_epi32(8);
   int u= 1;
   for(; u + 7 <= n; u+= 8, id= _mm256_add_epi32(id, eight)) {
    const __m256i w= ld(wr + u);
    const __m256i act= _mm256_and_si256(_mm256_cmpgt_epi32(w, zero), _mm256_cmpeq_epi32(ld(sv.data() + u), m1));
    const __m256i val= _mm256_sub_epi32(_mm256_add_epi32(ld(lb.data() + u), vb), _mm256_add_epi32(w, w));
    const __m256i better= _mm256_and_si256(act, _mm256_cmpgt_epi32(bv, val));
    bv= _mm256_blendv_epi8(bv, val, better), bi= _mm256_blendv_epi8(bi, id, better);
   }
   V b= INF;
   int bx= 0;
   reduce_min(bv, bi, b, bx);
   for(; u <= n; ++u)
    if(wr[u] > 0 && sv[u] == -1)
     if(const V val= lb[u] + base - 2 * wr[u]; val < b) b= val, bx= u;
   if(bx) key[x]= b, slack[x]= bx, kst[x]= clk;
  } else
   for(int u= 1; u <= n; ++u)
    if(const E& ex= row(x, u); w0(ex.u, ex.v) > 0 && st[u] != x && S[st[u]] == 0) update_slack(u, x);
 }
 V slack_value(int x) const { return key[x] - cf(x) * T; }
 void q_push(int x) {
  if(x <= n) q.push_back(x);
  else
   for(int y: flo[x]) q_push(y);
 }
 void set_st(int x, int b) {
  st[x]= b;
  if(x > n)
   for(int y: flo[x]) set_st(y, b);
 }
 int get_pr(int b, int xr) {
  int pr= std::find(flo[b].begin(), flo[b].end(), xr) - flo[b].begin();
  if(pr % 2 == 1) {
   std::reverse(flo[b].begin() + 1, flo[b].end());
   return (int)flo[b].size() - pr;
  }
  return pr;
 }
 void set_match(int u, int v) {
  const E r= rep(u, v);
  match[u]= r.v;
  if(u <= n) return;
  const int xr= child_of(u, r.u), pr= get_pr(u, xr);
  for(int i= 0; i < pr; ++i) set_match(flo[u][i], flo[u][i ^ 1]);
  set_match(xr, v);
  std::rotate(flo[u].begin(), flo[u].begin() + pr, flo[u].end());
 }
 void augment(int u, int v) {
  for(;;) {
   const int xnv= st[match[u]];
   set_match(u, v);
   if(!xnv) return;
   set_match(xnv, st[pa[xnv]]);
   u= st[pa[xnv]], v= xnv;
  }
 }
 int get_lca(int u, int v) {
  for(++tvis; u || v; std::swap(u, v)) {
   if(u == 0) continue;
   if(vis[u] == tvis) return u;
   vis[u]= tvis;
   u= st[match[u]];
   if(u) u= st[pa[u]];
  }
  return 0;
 }
 void add_blossom(int u, int lca, int v) {
  int b= n + 1;
  while(b <= nx && st[b]) ++b;
  if(b > nx) ++nx;
  match[b]= match[lca];
  flo[b].clear();
  flo[b].push_back(lca);
  for(int x= u, y; x != lca; x= st[pa[y]]) flo[b].push_back(x), flo[b].push_back(y= st[match[x]]), q_push(y);
  std::reverse(flo[b].begin() + 1, flo[b].end());
  for(int x= v, y; x != lca; x= st[pa[y]]) flo[b].push_back(x), flo[b].push_back(y= st[match[x]]), q_push(y);
  // 中に入る花の双対は止め、中の頂点はすべて偶にする。新しい花の双対は 0 から始める。
  for(int x: flo[b]) set_blossom_coef(x, 0), set_vertex_coef(x, -1), S[x]= -1, bp[x]= b;
  rt[b]= rt[lca], S[b]= 0, sv[b]= 2, lb[b]= -2 * T, bp[b]= 0;
  set_st(b, b);
  // 花の行。どの子の頂点も偶なので、元の頂点 x への辺 (m, x) の縮約費用は lb[m] - T + lab(x) - 2w(m, x)。
  E* rb= &row(b, 0);
  std::fill(rb, rb + N2, E{0, 0});
  std::fill(bd.begin(), bd.begin() + nx + 1, INF);
  const __m256i vT= _mm256_set1_epi32(T), zero= _mm256_setzero_si256(), one32= _mm256_set1_epi32(1), m1= _mm256_set1_epi32(-1);
  const __m256i vn1= _mm256_set1_epi32(n1), seq= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
  for(int xs: flo[b]) {
   const E* rx= xs > n ? &row(xs, 0) : nullptr;
   const __m256i vxs= _mm256_set1_epi32(xs), lbxs= _mm256_set1_epi32(xs <= n ? lb[xs] : 0);
   int x= 1;
   for(; x + 7 <= n; x+= 8) {
    const __m256i id= _mm256_add_epi32(_mm256_set1_epi32(x), seq);
    __m256i m, w, lm;
    if(!rx) m= vxs, w= ld(wv.data() + (size_t)xs * n1 + x), lm= lbxs;
    else {
     // 花 xs の行の 8 本の辺の、花の中の端 (E の u) を集める。
     m= _mm256_i32gather_epi32((const int*)(rx + x), _mm256_add_epi32(seq, seq), 4);
     w= _mm256_i32gather_epi32(wv.data(), _mm256_add_epi32(_mm256_mullo_epi32(m, vn1), id), 4);
     lm= _mm256_i32gather_epi32(lb.data(), m, 4);
    }
    const __m256i s4= ld(sv.data() + x);
    const __m256i labx= _mm256_sub_epi32(_mm256_add_epi32(ld(lb.data() + x), _mm256_and_si256(_mm256_cmpeq_epi32(s4, one32), vT)), _mm256_and_si256(_mm256_cmpeq_epi32(s4, m1), vT));
    const __m256i d= _mm256_sub_epi32(_mm256_add_epi32(_mm256_sub_epi32(lm, vT), labx), _mm256_add_epi32(w, w));
    const __m256i cur= ld(bd.data() + x);
    const __m256i better= _mm256_and_si256(_mm256_cmpgt_epi32(w, zero), _mm256_cmpgt_epi32(cur, d));
    if(_mm256_testz_si256(better, better)) continue;
    _mm256_storeu_si256((__m256i*)(bd.data() + x), _mm256_blendv_epi8(cur, d, better));
    for(int mk= _mm256_movemask_ps(_mm256_castsi256_ps(better)); mk; mk&= mk - 1) {
     const int j= __builtin_ctz(mk);
     rb[x + j]= E{rx ? rx[x + j].u : xs, x + j};
    }
   }
   for(; x <= nx; ++x) {
    if(x == b) continue;
    const E c= rep(xs, x);
    if(!c.u || w0(c.u, c.v) == 0) continue;
    if(const V d= dlt(c); d < bd[x]) bd[x]= d, rb[x]= c;
   }
  }
  ftime[b]= ++fclk;
  set_slack(b);
 }
 void expand_blossom(int b) {
  for(int y: flo[b]) set_st(y, y);
  const int xr= child_of(b, rep(b, pa[b]).u), pr= get_pr(b, xr);
  for(int y: flo[b]) bp[y]= 0;
  for(int i= 0; i < pr; i+= 2) {
   const int xs= flo[b][i], xns= flo[b][i + 1];
   pa[xs]= rep(xns, xs).u;
   set_label(xs, 1), set_label(xns, 0), rt[xs]= rt[xns]= rt[b];
   reset_slack(xs), set_slack(xns);
   q_push(xns);
  }
  set_label(xr, 1), rt[xr]= rt[b], pa[xr]= pa[b];
  reset_slack(xr);
  for(size_t i= pr + 1; i < flo[b].size(); ++i) {
   const int xs= flo[b][i];
   set_label(xs, -1), set_slack(xs);
  }
  st[b]= 0, S[b]= -1;
 }
 // 増加路を流したら、両端の木の根を dead に入れて true を返す。
 bool on_found_edge(const E& x) {
  const int u= st[x.u], v= st[x.v];
  if(S[v] == -1) {
   pa[v]= x.u, set_label(v, 1), rt[v]= rt[u];
   const int nu= st[match[v]];
   reset_slack(v), reset_slack(nu);
   set_label(nu, 0), rt[nu]= rt[u], q_push(nu);
  } else if(S[v] == 0) {
   const int lca= get_lca(u, v);
   if(!lca) {
    dead.push_back(rt[u]), dead.push_back(rt[v]);
    return augment(u, v), augment(v, u), true;
   }
   add_blossom(u, lca, v);
  }
  return false;
 }
 // 流した 2 本の木を壊す。印を消してから、生きている木の偶の頂点からの slack を付け直し、きつい辺があれば伸ばす。
 void dissolve() {
  const int ra= dead[0], rb= dead[1];
  dead.clear();
  cand.clear();
  for(int x= 1; x <= nx; ++x)
   if(st[x] == x && S[x] != -1 && (rt[x] == ra || rt[x] == rb)) set_label(x, -1), reset_slack(x), cand.push_back(x);
  for(int x: cand) {
   if(S[x] != -1) continue;  // 先の頂点から伸びた木に拾われた
   set_slack(x);
   if(slack[x] && slack_value(x) == 0) on_found_edge(rep(slack[x], x));
  }
 }
 // 偶の頂点 u の行を走査する。増加路を流したら true。
 bool scan(int u) {
  const V* wr= &wv[u * n1];
  const V lu= lab(u);
  auto one= [&](int v) {
   if(wr[v] > 0 && st[v] == v && st[u] != v && S[v] != 1) {
    if(lu + lab(v) - 2 * wr[v] == 0) return on_found_edge(E{u, v});
    update_slack(u, v);
   }
   return false;
  };
  const __m256i vlu= _mm256_set1_epi32(lu), vT= _mm256_set1_epi32(T), zero= _mm256_setzero_si256();
  const __m256i one32= _mm256_set1_epi32(1), m1= _mm256_set1_epi32(-1);
  int v= 1;
  for(; v + 7 <= n; v+= 8) {
   const __m256i w= ld(wr + v), s4= ld(sv.data() + v), S4= ld(S.data() + v);
   const __m256i id= _mm256_add_epi32(_mm256_set1_epi32(v), _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7));
   const __m256i labv= _mm256_sub_epi32(_mm256_add_epi32(ld(lb.data() + v), _mm256_and_si256(_mm256_cmpeq_epi32(s4, one32), vT)), _mm256_and_si256(_mm256_cmpeq_epi32(s4, m1), vT));
   const __m256i dd= _mm256_sub_epi32(_mm256_add_epi32(vlu, labv), _mm256_add_epi32(w, w));
   const __m256i act= _mm256_andnot_si256(_mm256_cmpeq_epi32(S4, one32), _mm256_and_si256(_mm256_cmpgt_epi32(w, zero), _mm256_cmpeq_epi32(ld(st.data() + v), id)));
   if(!_mm256_testz_si256(act, _mm256_cmpeq_epi32(dd, zero))) {
    for(int j= 0; j < 8; ++j)
     if(one(v + j)) return true;
    continue;
   }
   const __m256i cand= _mm256_add_epi32(_mm256_add_epi32(dd, vT), _mm256_and_si256(_mm256_cmpeq_epi32(S4, zero), vT));
   const __m256i kv= ld(key.data() + v);
   const __m256i upd= _mm256_and_si256(act, _mm256_cmpgt_epi32(kv, cand));
   if(_mm256_testz_si256(upd, upd)) continue;
   _mm256_storeu_si256((__m256i*)(key.data() + v), _mm256_blendv_epi8(kv, cand, upd));
   for(int m= _mm256_movemask_ps(_mm256_castsi256_ps(upd)); m; m&= m - 1) slack[v + __builtin_ctz(m)]= u, kst[v + __builtin_ctz(m)]= clk;
  }
  for(; v <= n; ++v)
   if(one(v)) return true;
  for(int b= n + 1; b <= nx; ++b)
   if(st[b] == b && S[b] != 1 && st[u] != b && row(b, u).u) {
    if(dlt(row(b, u)) == 0) {
     if(on_found_edge(rep(u, b))) return true;
    } else update_slack(u, b);
   }
  return false;
 }
 std::vector<int> solve() {
  nx= n;
  for(int u= 0; u <= n; ++u) st[u]= u, flo[u].clear();
  V wmax= 0;
  for(int u= 1; u <= n; ++u)
   for(int v= 1; v <= n; ++v) wmax= std::max(wmax, w0(u, v));
  for(int u= 1; u <= n; ++u) lb[u]= wmax, rt[u]= u, pa[u]= 0, set_label(u, 0), q.push_back(u);
  for(;;) {
   while(qh < (int)q.size()) {
    const int u= q[qh++];
    if(S[st[u]] != 0) continue;
    if(scan(u)) dissolve();
   }
   // 双対の動かし方。slack は 2 倍の尺度で比べる (印の無い頂点への辺は縮約費用の 2 倍、偶どうしは縮約費用、奇の花は双対)。
   V best= INF;
   // 元の頂点の最小。最小を与えたものが古ければ付け直してやり直す (古い slack は本当の値以下)。
   for(;;) {
    const __m256i vT= _mm256_set1_epi32(T), zero= _mm256_setzero_si256(), one32= _mm256_set1_epi32(1), m1= _mm256_set1_epi32(-1);
    const __m256i vinf= _mm256_set1_epi32(INF), eight= _mm256_set1_epi32(8);
    __m256i bv= vinf, bi= zero, id= _mm256_setr_epi32(1, 2, 3, 4, 5, 6, 7, 8);
    int x= 1;
    for(; x + 7 <= n; x+= 8, id= _mm256_add_epi32(id, eight)) {
     const __m256i S4= ld(S.data() + x), kv= ld(key.data() + x);
     const __m256i act= _mm256_andnot_si256(_mm256_cmpeq_epi32(S4, one32), _mm256_and_si256(_mm256_cmpeq_epi32(ld(st.data() + x), id), _mm256_cmpgt_epi32(vinf, kv)));
     const __m256i val= _mm256_sub_epi32(_mm256_sub_epi32(kv, vT), _mm256_and_si256(_mm256_cmpeq_epi32(S4, zero), vT));
     const __m256i sc= _mm256_add_epi32(val, _mm256_and_si256(_mm256_cmpeq_epi32(S4, m1), val));
     const __m256i better= _mm256_and_si256(act, _mm256_cmpgt_epi32(bv, sc));
     bv= _mm256_blendv_epi8(bv, sc, better), bi= _mm256_blendv_epi8(bi, id, better);
    }
    V b= INF;
    int bx= 0;
    reduce_min(bv, bi, b, bx);
    for(; x <= n; ++x)
     if(st[x] == x && S[x] != 1 && slack[x])
      if(const V val= slack_value(x) * (S[x] == -1 ? 2 : 1); val < b) b= val, bx= x;
    if(bx && !slack_valid(bx)) {
     set_slack(bx);
     continue;
    }
    best= b;
    break;
   }
   for(int x= n + 1; x <= nx; ++x) {
    if(st[x] != x || S[x] == 1 || !slack[x]) continue;
    V val= slack_value(x) * (S[x] == -1 ? 2 : 1);
    if(val >= best) continue;
    if(!slack_valid(x)) {
     set_slack(x);
     if(!slack[x]) continue;
     val= slack_value(x) * (S[x] == -1 ? 2 : 1);
     if(val >= best) continue;
    }
    best= val;
   }
   for(int b= n + 1; b <= nx; ++b)
    if(st[b] == b && S[b] == 1) best= std::min(best, lab(b));
   const V d= best / 2;
   // 偶の頂点の双対の最小 (lb - T)。これが d 以下なら終わり。偶の頂点が無ければ空いている頂点も無い。
   {
    const __m256i m1= _mm256_set1_epi32(-1);
    __m256i mv= _mm256_set1_epi32(INF);
    int u= 1;
    for(; u + 7 <= n; u+= 8) {
     const __m256i l= ld(lb.data() + u);
     const __m256i act= _mm256_and_si256(_mm256_cmpeq_epi32(ld(sv.data() + u), m1), _mm256_cmpgt_epi32(mv, l));
     mv= _mm256_blendv_epi8(mv, l, act);
    }
    V mlb= INF;
    int dummy= 0;
    reduce_min(mv, _mm256_setzero_si256(), mlb, dummy);
    for(; u <= n; ++u)
     if(sv[u] == -1) mlb= std::min(mlb, lb[u]);
    if(mlb >= INF || best >= INF || mlb - T <= d) break;
   }
   T+= d;
   {
    const __m256i vT= _mm256_set1_epi32(T), zero= _mm256_setzero_si256(), one32= _mm256_set1_epi32(1);
    int x= 1;
    auto one= [&](int y) {
     if(st[y] == y && S[y] != 1 && slack[y] && slack_value(y) == 0 && slack_valid(y))
      if(on_found_edge(rep(slack[y], y))) dissolve();
    };
    for(; x + 7 <= n; x+= 8) {
     const __m256i S4= ld(S.data() + x);
     const __m256i val= _mm256_sub_epi32(_mm256_sub_epi32(ld(key.data() + x), vT), _mm256_and_si256(_mm256_cmpeq_epi32(S4, zero), vT));
     const __m256i hit= _mm256_andnot_si256(_mm256_cmpeq_epi32(S4, one32), _mm256_cmpeq_epi32(val, zero));
     for(int m= _mm256_movemask_ps(_mm256_castsi256_ps(hit)); m; m&= m - 1) one(x + __builtin_ctz(m));
    }
    for(; x <= nx; ++x) one(x);
   }
   for(int b= n + 1; b <= nx; ++b)
    if(st[b] == b && S[b] == 1 && lab(b) == 0) expand_blossom(b);
  }
  std::vector<int> mt(n, -1);
  for(int u= 1; u <= n; ++u)
   if(match[u]) mt[u - 1]= match[u] - 1;
  return mt;
 }
};
using i64= long long;
// 重みが正で自己ループの無い辺 {u, v, w} の、重みの和が最大のマッチングの相手の配列。
inline std::vector<int> solve(int n, const std::vector<std::tuple<int, int, i64>>& e) {
 i64 wmax= 0;
 for(auto& [u, v, w]: e) wmax= std::max(wmax, w);
 if(wmax * 8 < (i64(1) << 30)) {
  Solver32 s(n);
  for(auto& [u, v, w]: e) s.add_edge(u, v, (int)w);
  return s.solve();
 }
 Solver64 s(n);
 for(auto& [u, v, w]: e) s.add_edge(u, v, w);
 return s.solve();
}
}  // namespace weighted_matching_internal
template <class W> WeightedMatching<W> max_weight_matching(int n, const std::vector<std::tuple<int, int, W>>& edges) {
 static_assert(std::is_integral_v<W>, "重みは整数の型");
 using i64= long long;
 std::vector<std::tuple<int, int, i64>> e;
 e.reserve(edges.size());
 for(auto& [u, v, w]: edges) {
  assert(0 <= u && u < n && 0 <= v && v < n);
  if(u != v && w > 0) e.emplace_back(u, v, (i64)w);
 }
 WeightedMatching<W> ret{W(0), weighted_matching_internal::solve(n, e)};
 // 組ませた 2 頂点の間の辺のうち、いちばん重いもの (核は多重辺の重い方を使う)。
 std::vector<i64> best(n, 0);
 for(auto& [u, v, w]: e)
  if(ret.mate[u] == v) best[u]= best[v]= std::max(best[u], w);
 for(int u= 0; u < n; ++u)
  if(ret.mate[u] > u) ret.weight+= (W)best[u];
 return ret;
}
template <class W> std::optional<WeightedMatching<W>> min_weight_perfect_matching(int n, const std::vector<std::tuple<int, int, W>>& edges) {
 static_assert(std::is_integral_v<W>, "重みは整数の型");
 using i64= long long;
 if(n & 1) return std::nullopt;
 if(n == 0) return WeightedMatching<W>{W(0), {}};
 i64 lo= std::numeric_limits<i64>::max(), hi= std::numeric_limits<i64>::min();
 for(auto& [u, v, w]: edges) {
  assert(0 <= u && u < n && 0 <= v && v < n);
  if(u != v) lo= std::min(lo, (i64)w), hi= std::max(hi, (i64)w);
 }
 if(lo > hi) return std::nullopt;  // 自己ループのほかに辺が無い
 const i64 C= (i64)(n / 2) * (hi - lo) + lo + 1;
 std::vector<std::tuple<int, int, i64>> e;
 e.reserve(edges.size());
 for(auto& [u, v, w]: edges)
  if(u != v) e.emplace_back(u, v, C - (i64)w);
 std::vector<int> mate= weighted_matching_internal::solve(n, e);
 for(int v= 0; v < n; ++v)
  if(mate[v] < 0) return std::nullopt;
 // 組ませた 2 頂点の間の辺のうち、いちばん軽いもの。
 std::vector<i64> best(n, std::numeric_limits<i64>::max());
 for(auto& [u, v, w]: edges)
  if(u != v && mate[u] == v) best[u]= best[v]= std::min(best[u], (i64)w);
 W sum= 0;
 for(int u= 0; u < n; ++u)
  if(mate[u] > u) sum+= (W)best[u];
 return WeightedMatching<W>{sum, std::move(mate)};
}
