/** \file OMPSimpleVector.h*/


#ifndef OMP_SIMPLEVECTOR_H_
#define OMP_SIMPLEVECTOR_H_

typedef struct  OMPCACell OMPCACell;

//! A vector structure abstraction.
/** The OMPSimpleVector holds integers and has a specific size. It has abstraction of pop and push (thus it works like a stack)
* and it supports atomic inserts.
*/
typedef struct OMPSimpleVector
{
	int m_size;
	int m_data[64];	
} OMPSimpleVector;

//! Insert an element to the vector.
/** Pushes an integer to the back of the vector if there is enough space.
* Returns index of position.
\param v The vector structure
\param element The element to insert
*/
int push_backsv(OMPSimpleVector* v, const int element) {
	int previousSize = v->m_size;
	(v->m_size)++;
	if(previousSize<64)	{
		v->m_data[previousSize] = element;
		return previousSize;
	} else {
		printf (":(\n");	
		--(v->m_size);
		return -1;
	}
}

//! Extract an element from the vector.
/** Pops an element from the back of the vector if any.
\param v The vector structure
\param cell The pointer to the item in which the extracted value is written
*/
int pop_backsv(OMPSimpleVector* v, int* cell) {
	if(v->m_size > 0) {
		int previousSize = (v->m_size)--;
		*cell = v->m_data[previousSize-1];
		return 1;
	} else {
		return -1;
	}
}

//! Insert an element to the vector in a thread-safe way.
/** Pushes an integer to the back of the vector if there is enough space. Uses atomic operations to make it thread-safe. 
* Returns index of position.
\param v The vector structure
\param element The element to insert
*/
int push_back_tssv(OMPSimpleVector* v, const int element) {
	int previousSize;
	
	#pragma omp atomic capture
	previousSize = (v->m_size)++;
	if(previousSize<48) {
		v->m_data[previousSize] = element;
		return previousSize;
	} else {
		#pragma omp atomic
		--(v->m_size);
		return -1;
	}
}

//! Set the vector to empty.
void resetsv(OMPSimpleVector* v) {
	v->m_size = 0;
}

//! Get number of elements in vector.
int sizesv(const OMPSimpleVector* v) {
	return v->m_size;
}

#endif
