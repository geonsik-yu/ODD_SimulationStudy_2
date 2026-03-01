// ===========================================================================
// <<  FUNCTION DEFINITIONS (PROGRAMMING UTILITY FUNCTIONS)  >>
// ===========================================================================

// <vector> <string> included in the header.
#include <algorithm>    // For std::sort, std::shuffle, std::count
#include <random>       // For std::mt19937, std::random_device 
#include <utility>      // For std::pair
#include <cstdlib>      // For rand(), RAND_MAX
#include "include/utils.h"

// 1) Joins multiple strings in a string vector:
void join(const std::vector<std::string>& v, char c, std::string& s) {
    s.clear();
    for (std::vector<std::string>::const_iterator p = v.begin(); p!= v.end(); ++p) 
        { s += *p;   if (p != v.end() - 1) s += c; }
}

// 2) Random Number Generation: returns random number generating function with Unif(0,1)
double rand01(std::default_random_engine& generator) {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}
// double rand01() { return ((double)rand() / (RAND_MAX)); } // Deprecated c-style function


// 3-1) Find agent indices with the given attribute vector OBJ & the directed quantile LB & UB.
std::vector<int> find_quantile_indices(const std::vector<double>& OBJ, int lb_quantile, int ub_quantile) {
    // -- Input vector should not be empty
    // -- Input "quantile" must be one of {1, 2, 3, 4}
    int n = static_cast<int> (OBJ.size());
    // Pair values with their original indices
    std::vector<std::pair<double, int>> value_index_pairs;
    for (int i = 0; i < n; ++i) {  value_index_pairs.emplace_back(OBJ[i], i); } // First -- OBJ value,  Second -- Index
    // Sort the pairs by value
    std::sort(value_index_pairs.begin(), value_index_pairs.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    // Determine range for the specified quantile
    int start = 0;       int end = 0;
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if      (lb_quantile == 0) {  start = 0;           } // Q0 (Min)
    else if (lb_quantile == 1) {  start = (n / 4);     } // Q1
    else if (lb_quantile == 2) {  start = (n / 2);     } // Q2
    else if (lb_quantile == 3) {  start = (3 * n / 4); } // Q3
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if      (ub_quantile == 1) {  end   = (n / 4);     } // Q1
    else if (ub_quantile == 2) {  end   = (n / 2);     } // Q2
    else if (ub_quantile == 3) {  end   = (3 * n / 4); } // Q3
    else if (ub_quantile == 4) {  end   = n;           } // Q4 (Max)
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Collect indices for the specified quantile
    std::vector<int> quantile_indices;
    for (int i = start; i < end; ++i) { quantile_indices.push_back(value_index_pairs[i].second);}
    return quantile_indices;
}

// 3-2) Find agent indices with the given attribute vector OBJ & the directed percentile chunk LB & UB.
std::vector<int> find_10_percentile_indices(const std::vector<double>& OBJ, int nth_chunk) {
    // -- Input vector should not be empty
    // -- Input "quantile" must be one of {0, 1, 2, 3, 4, 5, ..., 9} >> 10 chunks
    int n = static_cast<int> (OBJ.size());
    // Pair values with their original indices
    std::vector<std::pair<double, int>> value_index_pairs;
    for (int i = 0; i < n; ++i) {  value_index_pairs.emplace_back(OBJ[i], i); } // First -- OBJ value,  Second -- Index
    // Sort the pairs by value
    std::sort(value_index_pairs.begin(), value_index_pairs.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    // Determine range for the specified quantile
    int start = ( n * nth_chunk / 10 );       
    int end   = ( n * (nth_chunk + 1) / 10 );
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Collect indices for the specified quantile
    std::vector<int> result_indices;
    for (int i = start; i < end; ++i) { result_indices.push_back(value_index_pairs[i].second);}
    return result_indices;
}

// 4) Random selection function for initially infectious population selection.
//std::vector<int> randomSelectElements(const std::vector<int>& inputVector, int N, std::mt19937& generator) { // Alternative with Mersenne Twister
std::vector<int> randomSelectElements(const std::vector<int>& inputVector, int N, std::default_random_engine& generator) {    
    std::vector<int> result;
    // Create a copy of the input vector
    std::vector<int> copyVector = inputVector;
    // Shuffle the elements in the copy vector
    std::shuffle(copyVector.begin(), copyVector.end(), generator);
    // Select the first N elements from the shuffled copy vector
    for (int i = 0; i < N; ++i) { result.push_back(copyVector[i]); }
    return result;
}

// 5) Count the number of occurrences of the target character in a char vector.
int countOccurrences(const std::vector<char>& characters, char target) {
    return std::count(characters.begin(), characters.end(), target);
}
