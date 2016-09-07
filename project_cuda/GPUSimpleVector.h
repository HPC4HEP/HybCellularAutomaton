/** \file GPUSimpleVector.h*/


#ifndef GPU_SIMPLEVECTOR_H_
#define GPU_SIMPLEVECTOR_H_

#include <cuda.h>

//! A vector structure abstraction.
/** The GPUSimpleVector holds integers and has a specific size. It has abstraction of pop and push (thus it works like a stack)
* and it supports atomic inserts.
*/
template<int maxSize, class T>
struct GPUSimpleVector
{
	__device__
	int push_back(const T& element)
	{

		auto previousSize = m_size;
		m_size++;
		if(previousSize<maxSize)
		{
			m_data[previousSize] = element;
			return previousSize;
		}
		else
		{
			--m_size;
			return -1;
		}
	};

	__device__
	int push_back_ts(const T& element)
	{
		auto previousSize = atomicAdd(&m_size, 1);
		if(previousSize<maxSize)
		{
			m_data[previousSize] = element;
			return previousSize;
		}
		else
		{
			atomicSub(&m_size, 1);
			return -1;
		}
	};

	__device__
	T pop_back()
	{
		if(m_size > 0)
		{
			auto previousSize = m_size--;
			return m_data[previousSize-1];
		}
		else
		return T();
	};

	__host__   __device__
	void reset()
	{
		m_size = 0;
	}
	;

	__host__   __device__
	int size() const
	{
		return m_size;
	}
	;
	int m_size;

	T m_data[maxSize];

};

#endif
