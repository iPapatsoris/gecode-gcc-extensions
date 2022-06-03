#ifndef H_DYNAMIC_PARTITION
#define H_DYNAMIC_PARTITION 

#include <vector>

using namespace std;

class DynamicPartition {
	vector<unsigned int> *elements;
	vector<unsigned int> *elementToPos;
	vector<bool> splitPoint;

	public:
	
	DynamicPartition(unsigned int totalElements) {
		elements = new vector<unsigned int>();
		elements->reserve(totalElements);
		elementToPos = new vector<unsigned int>();
		elementToPos->reserve(totalElements);
		for (unsigned int i = 0; i < totalElements; i++) {
			elements->push_back(i);
			elementToPos->push_back(i);
			splitPoint.push_back(0);
		}
	}

	DynamicPartition() {}
};

#endif