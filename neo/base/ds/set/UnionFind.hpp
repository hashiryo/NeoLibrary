#pragma once
#include <utility>
#include <vector>
// 素集合の合併。根は p[x] < 0 で、その絶対値が大きさ。合併は大きいほうに寄せ、find で経路を縮める。
struct UnionFind {
 std::vector<int> p;
 explicit UnionFind(int n= 0): p(n, -1) {}
 int find(int x) { return p[x] < 0 ? x : p[x]= find(p[x]); }
 bool same(int a, int b) { return find(a) == find(b); }
 void unite(int a, int b) {
  a= find(a), b= find(b);
  if(a == b) return;
  if(p[a] > p[b]) std::swap(a, b);
  p[a]+= p[b], p[b]= a;
 }
 int size(int x) { return -p[find(x)]; }
};
