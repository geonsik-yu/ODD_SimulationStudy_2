#ifndef UTILS_H
#define UTILS_H

#include <vector>       
#include <string>       
#include <random>       

// 1) Joins multiple strings in a string vector:
void join(const std::vector<std::string>& v, char c, std::string& s);

// 2) Random Number Generation: returns random number generating function with Unif(0,1)
double rand01(std::default_random_engine& generator);

// 3-1) Find agent indices with the given attribute vector OBJ & the directed quantile LB & UB.
std::vector<int> find_quantile_indices(const std::vector<double>& OBJ, int lb_quantile, int ub_quantile);
// 3-2) Find agent indices with the given attribute vector OBJ & the directed percentile chunk LB & UB.
std::vector<int> find_10_percentile_indices(const std::vector<double>& OBJ, int nth_chunk);

// 4) Random selection function for initially infectious population selection.
std::vector<int> randomSelectElements(const std::vector<int>& inputVector, int N, std::default_random_engine& generator);   

// 5) Count the number of occurrences of the target character in a char vector.
int countOccurrences(const std::vector<char>& characters, char target);

#endif
