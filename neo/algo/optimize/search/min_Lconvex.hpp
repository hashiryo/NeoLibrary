#pragma once
#include <utility>
#include <vector>
// L 凸な f の最小を、座標を alpha ずつ動かす山登りで求める。alpha は改善しなくなるまで半分に
// していく。返すのは (f の最小値, その x)。座標の数は高々 30 程度を想定する (2^n の部分集合を試す)。
// 今の Library の min_Lconvex を写したもの。
template <class TD, class TR, class F> std::pair<TR, std::vector<TD>> min_Lconvex(const F& f, std::vector<TD> x, TD alpha) {
 TR f0= f(x), f1= f0, fS;
 for(int n= x.size(); alpha; f0 == f1 ? alpha>>= 1 : f0= f1) {
  std::vector<TD> x0{x};
  for(int S= 1; S < (1 << n) - 1; S++) {
   std::vector<TD> xS{x0};
   for(int i= 0; i < n; i++)
    if((S >> i) & 1) xS[i]+= alpha;
   if((fS= f(xS)) < f1) f1= fS, x= std::move(xS);
  }
 }
 return {f1, std::move(x)};
}
