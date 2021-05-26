#ifndef VEC_FUNCS_H
#define VEC_FUNCS_H

#include <numeric>

namespace {

template <typename T>
double getVectorMean(const std::vector<T> &vec) {
	return std::accumulate(vec.begin(), vec.end(), 0.0)/vec.size();
};

template <typename T> 
T getVectorMax(const std::vector<T> &vec) {
	return *std::max_element(vec.begin(), vec.end());
};

template <typename T> 
T getVectorMin(const std::vector<T> &vec) {
	return *std::min_element(vec.begin(), vec.end());
};

template <typename T> 
int getVectorArgMax(const std::vector<T> &vec) {
	return std::distance(vec.begin(), std::max_element(vec.begin(), vec.end()));
};

template <typename T> 
int getVectorArgMin(const std::vector<T> &vec) {
	return std::distance(vec.begin(), std::min_element(vec.begin(), vec.end()));
};

template <typename T>
T getVectorMedian(const std::vector<T> &vec) {
	std::vector<T> vecCopy=vec;  
	std::nth_element(vecCopy.begin(), vecCopy.begin()+vecCopy.size()/2, vecCopy.end());
	T median=vecCopy[vecCopy.size()/2];
	return median;
};


template <typename T>
std::vector<T> linspace(T min, T max, int num) {
	std::vector<T> vec;
	for (int i=0; i < num; i++ ) {
		vec.push_back((max-min)/num*i);
	};
	return vec;
};

template <typename T>
std::vector<T> arange(int min,int max) {
	std::vector<T> vec;
	for (int i=min; i < max; i++ ) {
		vec.push_back(T(i));
	};
	return vec;
};

int indexToBytes(int index) {
	return index/3*2*4;
};



template<typename T>
std::vector<double> finiteDifference(const std::vector<T> &vec) {
	std::vector<double> derivative(vec.size());
	for (int i=1; i < vec.size(); i++) {
		derivative[i]=vec[i]-vec[i-1];
	}
	derivative[0]=derivative[1];
	return derivative;
};


template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &orig)
{   
	std::vector<T> ret;
	for(auto v=orig.begin(); v!=orig.end(); v++) {
		ret.insert(ret.end(), v->begin(), v->end());                                                                                   }
	return ret;
}   

int roundUp(int numToRound, int multiple)
{
	if (multiple == 0)
		return numToRound;

	int remainder = abs(numToRound) % multiple;
	if (remainder == 0)
		return numToRound;

	if (numToRound < 0)
		return -(abs(numToRound) - remainder);
	else
		return numToRound + multiple - remainder;
}
}

#endif
