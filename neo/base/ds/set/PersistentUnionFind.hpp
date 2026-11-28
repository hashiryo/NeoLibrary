#pragma once
#include <utility>
#include "neo/base/ds/range/PersistentArray.hpp"
// 永続 union-find。親の配列を永続配列 (M = 64) に持ち、unite は新しい版を返す。経路圧縮はしない。
// 合併は大きさの大きいほうに寄せるので、find は O(log n) 段、各段の get は O(log_64 n) ですむ。
// 版の写しは根 1 個。UnionFind と同じく、根は p[x] < 0 で、その絶対値が大きさ。
struct PersistentUnionFind {
 PersistentArray<int, 64> p;
 PersistentUnionFind()= default;
 explicit PersistentUnionFind(int n): p(n, -1) {}
 int find(int x) const {
  for(int y; (y= p.get(x)) >= 0;) x= y;
  return x;
 }
 bool same(int a, int b) const { return find(a) == find(b); }
 int size(int x) const { return -p.get(find(x)); }
 PersistentUnionFind unite(int a, int b) const {
  a= find(a), b= find(b);
  if(a == b) return *this;
  PersistentUnionFind res= *this;
  if(p.get(a) > p.get(b)) std::swap(a, b);
  res.p= res.p.set(b, a);
  res.p= res.p.set(a, p.get(a) + p.get(b));
  return res;
 }
 static void reserve(std::size_t n) { PersistentArray<int, 64>::reserve(n); }
};
