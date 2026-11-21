#pragma once
#include "neo/base/core/set_fps/generic.hpp"
#include "neo/base/core/set_fps/ranked_xor.hpp"
// power projection と、合成の転置。f の長さは 2 冪 (2^n)。
//   power_projection(f, w, M): ans_m = sum_S w[S] [S] f^m (m = 0..M-1)。f[∅] は何でもよい。
//   egf_T(f): [X^{全体}] f^k / k! (k = 0..n)。f[∅] は 0 とみなす。
//   egf_T(f, g): [X^{全体}] f^k / k! g (k = 0..n)。f[∅] は 0 とみなす。
//   c0 = f[∅]、g = f - c0、V_k = sum_S w[S] (g^k / k!)[S] とすると ans_m = sum_k (m)_k c0^{m-k} V_k ((m)_k は下降階乗) で、
//   V は合成 (composite.hpp) の転置を w に当てて求める。ans は u_k(m+1) = c0 u_k(m) + k u_{k-1}(m) で進めるので、c0 は何でもよい。
//   法が奇数の ZMod (5 <= n <= 22) は rank 付きの XOR 変換の速い核、ほかの環はランク付きのゼータ変換。V は O(n^2 2^n)、ans は O(n M)。
namespace sps {
namespace sps_pp_internal {
// V_k = sum_S w[S] (g^k / k!)[S]。g[∅] は読まない。
template <class T> std::vector<T> egf_proj(const std::vector<T>& g, const std::vector<T>& w) {
 const int n= sps_xor::log2_size<T>(g.size());
 if constexpr(sps_xor::is_zmod<T>::value) {
  if(sps_xor::fast<T>(n, 5)) {
   const sps_xor::MK k(T::mod());
   std::vector<sps_xor::u32> v;
   sps_xor::dispatch<5>(n, [&]<int N>() { v= sps_xor::egf_proj_n<N>(sps_xor::raw(g), sps_xor::raw(w), k); });
   std::vector<T> V(n + 1);
   for(int i= 0; i <= n; ++i) V[i]= T::raw(v[i]);
   return V;
  }
 }
 const std::size_t N= g.size();
 if(N == 1) return {w[0]};  // 今の Library の egf_T(f, g) は n = 0 で g[1] を読むので、ここで返す
 std::vector<T> f(g), wr(N);
 f[0]= T();
 for(std::size_t s= 0; s < N; ++s) wr[s]= w[(N - 1) ^ s];
 return sps_generic::egf_T(f, wr);
}
}  // namespace sps_pp_internal
template <class T> std::vector<T> egf_T(const std::vector<T>& f) {
 std::vector<T> w(f.size());
 w.back()= T(1);
 return sps_pp_internal::egf_proj(f, w);
}
template <class T> std::vector<T> egf_T(const std::vector<T>& f, const std::vector<T>& g) {
 const std::size_t N= f.size();
 assert(N == g.size());
 std::vector<T> w(N);
 for(std::size_t s= 0; s < N; ++s) w[s]= g[(N - 1) ^ s];
 return sps_pp_internal::egf_proj(f, w);
}
template <class T> std::vector<T> power_projection(const std::vector<T>& f, const std::vector<T>& w, int M) {
 assert(f.size() == w.size());
 const std::vector<T> V= sps_pp_internal::egf_proj(f, w);
 const int n= int(V.size()) - 1;
 const T c0= f[0];
 std::vector<T> ans(M);
 if constexpr(sps_xor::is_zmod<T>::value) {
  if((T::mod() & 1) && T::mod() >= 3 && n < 24) {
   std::vector<sps_xor::u32> v(n + 1);
   for(int i= 0; i <= n; ++i) v[i]= V[i].val();
   sps_xor::proj_post(v, c0.val(), M, sps_xor::raw(ans), sps_xor::MK(T::mod()));
   return ans;
  }
 }
 std::vector<T> u(n + 1);
 u[0]= T(1);
 for(int m= 0; m < M; ++m) {
  T s= T();
  for(int i= 0; i <= n; ++i) s+= V[i] * u[i];
  ans[m]= s;
  for(int i= n; i >= 1; --i) u[i]= u[i] * c0 + u[i - 1] * T(i);
  u[0]*= c0;
 }
 return ans;
}
}  // namespace sps
