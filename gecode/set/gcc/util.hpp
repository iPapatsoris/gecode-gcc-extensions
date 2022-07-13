#ifndef H_UTIL
#define H_UTIL

#include <unordered_map>
#include <unordered_set>

using namespace std;

typedef unordered_map<int, unordered_set<int>> MapToSet;

struct EdgeInfo {
	int src;
	int dest;
	bool lowerBoundViolation;
	int violationBound;

	EdgeInfo(int src, int dest, bool lowerBoundViolation = false, int violationBound = -1) 
		: src(src), dest(dest), lowerBoundViolation(lowerBoundViolation),
			violationBound(violationBound) {}
	EdgeInfo() {}
};

#endif