#ifndef H_UTIL
#define H_UTIL

#include <unordered_map>
#include <unordered_set>

using namespace std;

typedef unordered_map<int, unordered_set<int>> MapToSet;

struct EdgeInfo {
	int src;
	int dest;

	EdgeInfo(int src, int dest) : src(src), dest(dest) {}
	EdgeInfo() {}
};

#endif