/*
 * BDDTricks.hh
 *
 *  created on: 03.05.2017
 *      author: M.Khaled
 */

#ifndef VECTORTRICKS_HH_
#define VECTORTRICKS_HH_

#include <array>
#include <vector>

#include "VectorTricks.hh"
#include "cuddObj.hh"

class VectorTricks {
public:

	/* constructs a one vector from a set of vectors by appending */
	template<typename T>
	static
	inline
	std::vector<T> UnrollVectors(const std::vector<std::vector<T>>& a){
		std::vector<T> v;
		for(size_t i=0; i<a.size(); i++)
			AppendVector(v,a[i]);
		return v;
	}

	/* removes the values of one vector from other  one*/
	template<typename T>
	static	
	inline
	void SubtractVector(std::vector<T>& a, const std::vector<T>& b){
		for(size_t i=0; i<b.size(); i++)
			a.erase(remove(a.begin(), a.end(), b[i]), a.end());
	}
	
	/* Appends vector's elements to another vector*/
	template<typename T>
	static	
	inline
	void AppendVector(std::vector<T>& a, const std::vector<T>& b){
		a.reserve(a.size() + b.size());
		a.insert(a.end(), b.begin(), b.end());
	}
	
	template<typename T>
	static	
	inline
	void AppendVectors(std::vector<T>& a, const std::vector<std::vector<T>>& b){
		for(size_t i=0; i<b.size(); i++)
			AppendVector(a, b[i]);
	}
	

	/* removes the elements in vector (a) not exsiting in the vector (b) */
	template<typename T>
	static
	inline
	void IntersectVector(std::vector<T>& a, const std::vector<T>& b){
		std::vector<size_t> marked4Removal;
		for(size_t i=0; i<a.size(); i++)
			if(!isVectorElement(b, a[i]))
				marked4Removal.push_back(i);

		for(int i=marked4Removal.size()-1; i>=0; i--)
			a.erase(a.begin()+marked4Removal[i]);
	}

	template<typename T>
	static	
	inline
	std::vector<T> ArrayToVector(const T* A, size_t size){
		std::vector<T> v(size);
		for (size_t i=0; i<size; i++)
			v[i]=A[i];
		return v;
	}
	
	/* Prints a vector to std-out or other ostream*/
	template<typename T>
	static	
	inline
	void PrintVector(const std::vector<T>& v, char separator = ',', bool do_endl=true, std::ostream& strm=std::cout){
		size_t size = v.size();
		for(size_t i=0; i<size; i++){
			strm << v[i];
			if(i < size-1 && separator != 'X') strm << separator;
		}
		if(do_endl)
			strm << std::endl;
	}
	
	template<typename T>
	static	
	inline
	void PrintArray(const T* v, size_t size, char separator = ',', bool do_endl=true, std::ostream& strm=std::cout){
		for (size_t i=0; i<size; i++)
		{
			strm << v[i];
			if(i<size-1 && separator != 'X') strm << separator;
		}
		if(do_endl)
			strm << std::endl;
	}
	
	/* checks whether an element is in the vector*/
	template<typename T>
	static	
	inline
	bool isVectorElement(const std::vector<T>& a, const T& e){
		if ((find(a.begin(),a.end(),e) != a.end()))
			return true;
		else
			return false;
	}
	
	/* checks whether a vector is included in another vector*/
	template<typename T>
	static	
	inline
	bool isVectorIncluded(const std::vector<T>& a, const std::vector<T>& b){
		for(size_t i=0; i<a.size(); i++)
			if (!isVectorElement(b,a[i]))
				return false;
	
		return true;
	}
	
	/* checks whether the vector intersects with another vector*/
	template<typename T>
	static	
	inline
	bool isVectorIntersects(const std::vector<T>& a, const std::vector<T>& b){
		for(size_t i=0; i<a.size(); i++)
			if (isVectorElement(b,a[i]))
				return true;
	
		return false;
	}
	
	/* checks whether a vector has has same elements of another vector*/
	template<typename T>
	static	
	inline
	bool isVectorEquivilant(const std::vector<T>& a, const std::vector<T>& b){
	
		if(isVectorIncluded(a,b) && isVectorIncluded(b,a))
			return true;
		else
			return false;
	}


};
#endif
