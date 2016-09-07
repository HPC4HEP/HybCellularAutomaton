/** \file OMPResultVector.h*/


#ifndef OMP_RESULTVECTOR_H_
#define OMP_RESULTVECTOR_H_

typedef struct int4 {
	int elem[4];
} int4;

//! A vector structure abstraction.
/** The OMPResultVector holds integer quadruplets and has a specific size. It has abstraction of pop and push (thus it works like a stack)
* and it supports atomic inserts.
*/
typedef struct OMPResultVector
{
	int m_size;
	int4 m_data[1000];	
} OMPResultVector;


//! Insert an element to the vector.
/** Pushes a quadruplet to the back of the vector if there is enough space.
* Returns index of position.
\param v The vector structure
\param element The element to insert
*/
int push_backrv(OMPResultVector* v, const int4* element) {
	int previousSize = v->m_size;
	(v->m_size)++;
	if(previousSize<1000)	{
		v->m_data[previousSize] = *element;
		return previousSize;
	} else {
		--(v->m_size);
		return -1;
	}
}

//! Insert an element to the vector in a thread-safe way.
/** Pushes a quadruplet to the back of the vector if there is enough space. Uses atomic operations to make it thread-safe. 
* Returns index of position.
\param v The vector structure
\param element The element to insert
*/
int push_backtsrv(OMPResultVector* v, const int4* element) {
	int previousSize;
	#pragma omp atomic capture	
	previousSize = (v->m_size)++;
	if(previousSize<1000)	{
		v->m_data[previousSize] = *element;
		return previousSize;
	} else {
		#pragma omp atomic
		--(v->m_size);
		return -1;
	}
}

//! Extract an element from the vector.
/** Pops an element from the back of the vector if any.
\param v The vector structure
\param cell The pointer to the item in which the extracted value is written
*/
int pop_backrv(OMPResultVector* v, int4* cell) {
	if(v->m_size > 0) {
		int previousSize = (v->m_size)--;
		*cell = v->m_data[previousSize-1];
		return 1;
	} else {
		return -1;
	}
}

//! Set the vector to empty.
void resetrv(OMPResultVector* v) {
	v->m_size = 0;
}
//! Get number of elements in vector.
int sizerv(OMPResultVector* v) {
	return v->m_size;
}


#endif
