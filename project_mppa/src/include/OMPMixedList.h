/** \file OMPMixedList.h*/


#ifndef MIXEDLIST__H
#define MIXEDLIST__H


//! A structure that allocates memory as lists on an arena
/** The MixedList structure allocates pieces of memory by request and in
* an atomic manner by giving out memory indices. This way multiple lists can reside in the same buffer in a compact manner.
* Nodes are inserted at athe start of the list, and they used in a way similar to stacks.
*/
typedef struct MixedList {
	int size;
	int* buffer;
	int cap;
} MixedList;

//! Initialize a MixedList structure.
/** A buffer is allocated to the strcture as well as size limits.
\param list The structure to initialize
\param buf The buffer that is allocated to the structure
\param cap The maximum number of elements that can be stored
*/
void init_ml (MixedList* list, char* buf, int cap) {
	list->size = 0;
	list->cap = cap;
	list->buffer = (int*) buf;
}

//! Insert a new node in the MixedList strcture
/** This function allocates a node for inserting a value, links the node to a previous node and sets the value. Returns the identifier
* of the node.
\param list The structure to initialize
\param top The identifier of the previous node
\param val The value to be inserted
*/
int push_back_ml (MixedList* list, int top, int val) {
	int prevsize;
	#pragma omp atomic capture
	prevsize = (list->size)++;
	if (prevsize < list->cap) {
		list->buffer[2*prevsize] = top;
		list->buffer[2*prevsize+1] = val;
		return prevsize;
	} else {
		#pragma omp atomic
		(list->size)--;
		return -1;
	}
}

//! Insert a new node in the MixedList structure and update our exterior reference atomically.
/** This function allocates a node for inserting a value, links the node to a previous node and sets the value. Returns the identifier
* of the node. Is used for updating the reference to the start of the list in a safe manner.
\param list The structure to initialize
\param dest The pointer to the previous node's id. Will hold the new node's id afterwards
\param val The value to be inserted
*/
int push_back_mlts (MixedList* list, int* dest, int val) {
	int prevsize;
	#pragma omp atomic capture
	prevsize = (list->size)++;
	if (prevsize < list->cap) {
		#pragma omp atomic capture
		{list->buffer[2*prevsize] = *dest; *dest = prevsize;}
		list->buffer[2*prevsize+1] = val;
		return prevsize;
	} else {
		#pragma omp atomic
		(list->size)--;
		return -1;
	}
}

//! Get the successor of a node.
/** This function gets the identifier of a linked list node and retrieves the id of the successor
\param list The structure to initialize
\param top The identifier of the current node
*/
int next_ml (MixedList* list, int top) {
	return list->buffer[2*top];
}

//! Get the value of a node.
/** This function gets the identifier of a linked list node and retrieves the value associated with the node
\param list The structure to initialize
\param top The identifier of the current node
*/
int fetch_ml (MixedList* list, int top) {
	return list->buffer[2*top+1];
}


#endif












