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
	int buffer[12800];
} MixedList;

//! Initialize a MixedList structure.
/** The structure is set to empty.
\param list The structure to initialize

*/
void init_ml (MixedList* list) {
	list->size = 0;
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
	if (prevsize < 3200) {
		list->buffer[2*prevsize] = top;
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












