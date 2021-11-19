#ifndef H_SYM_UTIL
#define H_SYM_UTIL

#include <unordered_map>
#include <unordered_set>

template <class T1, class T2>
struct MapToSet {
	unordered_map<T1, unordered_set<T2>> map;
	MapToSet() {}
	MapToSet(const MapToSet &c) {
		map = c.map;
	}
};

#endif