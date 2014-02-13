#ifndef UTILS_VECTOR_UTILS_H_
#define UTILS_VECTOR_UTILS_H_

#include <vector>

using namespace std;

// Clear all memory allocated by this vector
template <typename T>
void ClearMemory(vector<T> & vt) {
    // Create an empty vector
    vector<T> emptyVector;
    // First clear the content
    vt.clear();
    // Then swap vt with the empty vector 
    vt.swap(emptyVector);
}

#endif
