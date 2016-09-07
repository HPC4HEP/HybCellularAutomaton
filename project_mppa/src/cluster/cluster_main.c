/*! \file cluster_main.c this file has the functions used by the compute cluster process. 
*This process retrieves problems from the i/o clusters and returns results
*/

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#include <HAL/hal/hal.h>
#include <mppaipc.h>
#include <mppa/osconfig.h>
#include <omp.h>

#include "shared_defs.h"
#include "settings.h"
#include "OMPResultVector.h"
#include "OMPHitsAndDoublets.h"
#include "OMPCACell.h"
#include "OMPSimpleVector.h"
#include "OMPMixedList.h"

//! The number of cores per cluster
#define TC 16

OMPResultVector results;

//! Retrieves quadruplets out of connected doublets
/*! This function gets the connected doublets and performs a DFS to produce the quadruplets.
\param doublets The data structure of the hit and doublet data
\param layer The current layer in the search
\param idx The identifier within the layer for the current doublet we are visiting
\param ml The lists of the connections between doublets
\param top An array that holds the indices of the start of every successor list
\param foundNtuplets The result returned for our problem
\param tmpNtuplet The current path of the search
\param minHitsPerNtuplet The number of layers involved
*/
void find_ntuplets(const OMPLayerDoublets* doublets, unsigned int layer, unsigned int idx, MixedList* ml, int** top, 
			OMPResultVector* foundNtuplets, OMPSimpleVector* tmpNtuplet, const unsigned int minHitsPerNtuplet) {
	int j;
	int4 found;
	int otherCell;

	if (layer == 0) {
		if (sizesv(tmpNtuplet) >= minHitsPerNtuplet - 1) {
			found.elem[0]=get_inner_hit_id(doublets, tmpNtuplet->m_data[2], 0);
			found.elem[1]=get_outer_hit_id(doublets, tmpNtuplet->m_data[2], 0);
			found.elem[2]=get_outer_hit_id(doublets, tmpNtuplet->m_data[1], 1);
			found.elem[3]=get_outer_hit_id(doublets, tmpNtuplet->m_data[0], 2);
			push_backtsrv(foundNtuplets, &found);
		}else
			return;
	} else {
		int ptr = top[layer][idx];
		while (ptr != -1) {
			int otherIdx = fetch_ml (ml, ptr);
			push_backsv(tmpNtuplet, otherIdx);
			find_ntuplets(doublets, layer-1, otherIdx, ml, top, foundNtuplets, tmpNtuplet, minHitsPerNtuplet);
			pop_backsv(tmpNtuplet, &otherCell);
			ptr = next_ml (ml, ptr);
		}
	}
}

//! everything the CA does is here.
/** The control flow is as follows:
* the process opens the file descriptors for communicating with IO, 
* loops so that it can handle the different request (didn't add a termination message) and executes the communication logic.
* It gets info on the size of the problem to receive (communication, also uses sync), allocates space for these layer pairs
* and starts reading. As soon as a layer pair is received, it performs the initialization of its cells. After we get more layer pairs,
* we connect the adjacent ones and finally we find the quadruplets. Then we push the results to IO.
*/
int
main(int argc __attribute__ ((unused)), char *argv[])
{
	int rank = 0;
	int i, j, status, idx1,idx2;
	float fargs[6];
	const int num_a = 6*sizeof(float)+3*(theNumberOfLayers-1)*sizeof(int);
	char args[num_a];
	int iargs[3*(theNumberOfLayers-1)];

	const char *root_sync = argv[1], *d_portal = argv[2], *a_portal = argv[3];
	const char *b_portal = argv[4], *c_portal = argv[5], *e_portal = argv[6];;
	

	/*Each cluster contributes a different bit to the root_sync mask.*/
	long long mask = (long long)1 << rank;

	/*Open the NoC special files.*/
	int root_sync_fd = mppa_open(root_sync, O_WRONLY);
	int d_portal_fd = mppa_open(d_portal, O_RDONLY);
	int a_portal_fd = mppa_open(a_portal, O_RDONLY);
	int b_portal_fd = mppa_open(b_portal, O_RDONLY);
	int c_portal_fd = mppa_open(c_portal, O_WRONLY);
	int e_portal_fd = mppa_open(e_portal, O_RDONLY);

	if (root_sync_fd < 0)
		printf ("Sync open error\n");
	if (a_portal_fd < 0)
                printf ("portal error\n");
	if (b_portal_fd < 0)
                printf ("portal error\n");
	if (c_portal_fd < 0)
                printf ("portal error\n");
	if (d_portal_fd < 0)
                printf ("portal error\n");
	if (e_portal_fd < 0)
                printf ("portal error\n");

	OMPLayerDoublets doublets[theNumberOfLayers-1];	
	

	int maxi = 0;
	int maxp = 0;
	/*arena to be managed*/
	char* buffer = malloc(1500000);
	
	int fail = 0;

	int num[3];
	/*initialize statically-sized asynchronous communications*/
	mppa_aiocb_t d_portal_aiocb[1] =
			{ MPPA_AIOCB_INITIALIZER(d_portal_fd, args, num_a) };
	mppa_aiocb_t c_portal_aiocb[1] =
			{ MPPA_AIOCB_INITIALIZER(c_portal_fd, &results, sizeof(OMPResultVector)) };

	mppa_aiocb_set_pwrite(c_portal_aiocb, &results, sizeof(OMPResultVector), 0);

	mppa_aiocb_set_trigger(d_portal_aiocb, 1);
	status |= mppa_aio_read(d_portal_aiocb);

	for (idx1 = 0; idx1 < 100; idx1++)
	for (idx2 = 0; idx2 < 7; idx2++) {

		char* origin = buffer;
		int left = 1500000;
		int flag = 0;
		unsigned int s;

		/*synchronize with io and gets parameters*/
		status |= mppa_write(root_sync_fd, &mask, sizeof(mask));
		if (idx1 != 0 || idx2 != 0)		
			mppa_aio_wait(c_portal_aiocb);
		status |= mppa_aio_wait(d_portal_aiocb);
		
		memcpy (fargs, args, 6*sizeof(float));
		memcpy (iargs, args+6*sizeof(float), 3*(theNumberOfLayers-1)*sizeof(int));

		num[0] = 2*sizeof(int)*iargs[0]+3*sizeof(float)*(iargs[1]+iargs[2]);
		num[1] = 2*sizeof(int)*iargs[3]+3*sizeof(float)*iargs[5];
		num[2] = 2*sizeof(int)*iargs[6]+3*sizeof(float)*iargs[8];

		//printf ("%d: Allocating memory\n", mppa_getpid());

		char* l0; char* l1; char* l2;

		int* liptr[3];
		
		/*allocate space*/
		left -= num[0];
		l0 = origin+left; /*allocated at the end of the buffer, is deallocated later*/	
				

		l1 = origin;
		left -= num[1];
		origin += num[1];

		l2 = origin;
		left -= num[2];
		origin += num[2];


		
		s = iargs[3];
		liptr[1] = (int*) origin;
		left -= s*sizeof(int);
		origin += s*sizeof(int);

		MixedList* ml = (MixedList*) origin;	
		left -=sizeof(MixedList);
		origin += sizeof(MixedList);
		if (left >= 0) {
			init_ml(ml, origin, 10000);
			left -= 2*10000*sizeof(int);
			origin += 2*10000*sizeof(int);
			if (left < 0)
				printf ("Out of memory: ml\n");
		} else {
			printf ("Out of memory: ml\n");
		}


		int* outerptr[2];
		MixedList* isOuterHitOfCell[2];

		isOuterHitOfCell[0] = (MixedList*) origin;
		left -= sizeof(MixedList);
		origin += sizeof(MixedList);
		if (left >= 0) {
			int size = iargs[0]+100;
			init_ml(isOuterHitOfCell[0], origin, size);
			left -= 2*size*sizeof(int);
			origin += 2*size*sizeof(int);
			if (left < 0)
				printf ("Out of memory: outer\n");
		} else {
			printf ("Out of memory: outer\n");
		}

		outerptr[0] = (int*) origin;
		left -= iargs[2]*sizeof(int);
		origin += iargs[2]*sizeof(int);
		if (left >= 0) {
			for (i = 0; i < iargs[2]; i++)
				outerptr[0][i] = -1;
		}

		isOuterHitOfCell[1] = (MixedList*) origin;
		left -= sizeof(MixedList);
		origin += sizeof(MixedList);
		if (left >= 0) {
			int size = iargs[3]+100;
			init_ml(isOuterHitOfCell[1], origin, size);
			left -= 2*size*sizeof(int);
			origin += 2*size*sizeof(int);
			if (left < 0)
				printf ("Out of memory: outer\n");
		} else {
			printf ("Out of memory: outer\n");
		}

		outerptr[1] = (int*) origin;
		left -= iargs[5]*sizeof(int);
		origin += iargs[5]*sizeof(int);
		if (left >= 0) {
			for (i = 0; i < iargs[5]; i++)
				outerptr[1][i] = -1;
		}

		int outer;
		int layer;
		int top, val;
	
		int thisCell[3];
		int otherCell[3];

		left -= 2*iargs[0]*sizeof(float);
		doublets[0].r = (float*) (origin+left);
		if (left < 0) {
			printf ("Out of memory: r1\n");
		}

		doublets[1].r = (float*) origin;
		left -= 2*iargs[3]*sizeof(float);	
		if (left < 0) {
			printf ("Out of memory: r2\n");
		}
		origin += 2*iargs[3]*sizeof(float);

			
		/*initialize hit and doublet receives
		get each layer in a different communication
		the idea is that we will wait for each layer just before we process it hiding other communications*/
		mppa_aiocb_t l0_portal_aiocb[1] =
			{ MPPA_AIOCB_INITIALIZER(a_portal_fd, l0, num[0]) };
		mppa_aiocb_t l1_portal_aiocb[1] =
			{ MPPA_AIOCB_INITIALIZER(b_portal_fd, l1, num[1]) };
		mppa_aiocb_t l2_portal_aiocb[1] =
			{ MPPA_AIOCB_INITIALIZER(e_portal_fd, l2, num[2]) };

		mppa_aiocb_set_trigger(l0_portal_aiocb, 1);
		status |= mppa_aio_read(l0_portal_aiocb);

		mppa_aiocb_set_trigger(l1_portal_aiocb, 1);
		status |= mppa_aio_read(l1_portal_aiocb);

		mppa_aiocb_set_trigger(l2_portal_aiocb, 1);
		status |= mppa_aio_read(l2_portal_aiocb);

		/*unlock the writes of io*/
		status |= mppa_write(root_sync_fd, &mask, sizeof(mask));

		status |= mppa_aio_wait(l0_portal_aiocb);

		doublets[0].size = iargs[0];
		doublets[0].indices = (int*) l0;
		l0 += 2*sizeof(int)*iargs[0];

		doublets[0].layers[0].size = iargs[1];
		doublets[0].layers[0].p = (float*) l0;
		l0 += 3*sizeof(float)*iargs[1];

		doublets[0].layers[1].size = iargs[2];
		doublets[0].layers[1].p = (float*) l0;

		/*Create layerpair 0-1*/
		if (left >= 0) {
			#pragma omp parallel for num_threads(TC) 
			for (j = 0; j < doublets[0].size; j++) {
				int in = doublets[0].indices[2*j];
				int out = doublets[0].indices[2*j+1];
				int inner = get_inner_hit_id (doublets, j, 0);
				float x = get_inner_x (doublets, inner, 0);
				float y = get_inner_y (doublets, inner, 0);
				doublets[0].r[2*j] = hypot(x,y);
				int outer = get_outer_hit_id (doublets, j, 0);
				x = get_outer_x (doublets, outer, 0);
				y = get_outer_y (doublets, outer, 0);
				doublets[0].r[2*j+1] = hypot(x,y);
				push_back_mlts (isOuterHitOfCell[0], &outerptr[0][outer], j);
			}
		}

		status |= mppa_aio_wait(l1_portal_aiocb);

		doublets[1].size = iargs[3];
		doublets[1].indices = (int*) l1;
		l1 += 2*sizeof(int)*iargs[3];

		doublets[1].layers[0].size = doublets[0].layers[1].size;
		doublets[1].layers[0].p = doublets[0].layers[1].p;
		doublets[1].layers[1].size = iargs[5];
		doublets[1].layers[1].p = (float*) l1;

		
		/*Create layerpair 1-2*/
		if (left >= 0) {
			#pragma omp parallel for num_threads(TC) 
			for (j = 0; j < doublets[1].size; j++) {
				int inner = get_inner_hit_id (doublets, j, 1);
				float x = get_inner_x (doublets, inner, 1);
				float y = get_inner_y (doublets, inner, 1);
				doublets[1].r[2*j] = hypot(x,y);
				int outer = get_outer_hit_id (doublets, j, 1);
				x = get_outer_x (doublets, outer, 1);
				y = get_outer_y (doublets, outer, 1);
				doublets[1].r[2*j+1] = hypot(x,y);
				push_back_mlts (isOuterHitOfCell[1], &outerptr[1][outer], j);
			}
	
			
		}
		/*connect 0-1-2*/
		if (left >= 0) {
			s = doublets[1].size;
			#pragma omp parallel for num_threads(TC) private(j, top, val, thisCell, otherCell)
			for (i = 0; i < s; i++) {
				top = -1;
				int inner = get_inner_hit_id(doublets, i, 1);
				thisCell[0] = i;
				thisCell[1] = inner;
				thisCell[2] = get_outer_hit_id(doublets, i, 1);
				
				int optr = outerptr[0][inner];
				/*loop through doublets sharing hit*/
				while (optr != -1) {
					otherCell[0] = fetch_ml (isOuterHitOfCell[0], optr);
					otherCell[1] = get_inner_hit_id(doublets, otherCell[0], 0);
					otherCell[2] = get_outer_hit_id(doublets, otherCell[0], 0);

					if (check_alignment_and_tag(doublets, thisCell, 1, otherCell,
								fargs[0], fargs[1], fargs[2],
								fargs[3], fargs[4], fargs[5])) {
						val = otherCell[0];
						top = push_back_ml (ml, top, val);
						if (top < 0) {
							printf ("Error: out of space\n");
							results.m_size = -1;
							left = -1;
							break;
						}
					}
					optr = next_ml (isOuterHitOfCell[0], optr);
				}
			
				liptr[1][i] = top;
			}
		}

		status |= mppa_aio_wait(l2_portal_aiocb);
		doublets[2].size = iargs[6];
		doublets[2].indices = (int*) l2;
		l2 += 2*sizeof(int)*iargs[6];

		doublets[2].layers[0].size = doublets[1].layers[1].size;
		doublets[2].layers[0].p = doublets[1].layers[1].p;

		doublets[2].layers[1].size = iargs[8];
		doublets[2].layers[1].p = (float*) l2;


		if (left < 0) {
			flag = 1;
			results.m_size = -1;
			goto res;
		}

		left += num[0];
		left += 2*iargs[0]*sizeof(float);

		

		doublets[2].r = (float*) origin;
		left -= 2*iargs[6]*sizeof(float);	
		origin += 2*iargs[6]*sizeof(float);

		s = iargs[6];
		liptr[2] = (int*) origin;
		left -= s*sizeof(int);
		origin += s*sizeof(int);

		if (left < 0) {
			printf ("Out of memory: r3\n");
			flag = 1;
			results.m_size = -1;
			goto res;
		}

		/*Create layerpair 2-3*/
		#pragma omp parallel for num_threads(TC) 
		for (j = 0; j < doublets[2].size; j++) {
			int inner = get_inner_hit_id (doublets, j, 2);
			float x = get_inner_x (doublets, inner, 2);
			float y = get_inner_y (doublets, inner, 2);
			doublets[2].r[2*j] = hypot(x,y);
			int outer = get_outer_hit_id (doublets, j, 2);
			x = get_outer_x (doublets, outer, 2);
			y = get_outer_y (doublets, outer, 2);
			doublets[2].r[2*j+1] = hypot(x,y);	
		}

		

		if (left < 0) {
			flag = 1;
			results.m_size = -1;
			goto res;
		}
		

		/*connect 1-2-3*/
		for (layer = 2; layer < theNumberOfLayers-1; layer++) {
			s = doublets[layer].size;
			#pragma omp parallel for num_threads(TC) private(j, top, val, thisCell, otherCell)
			for (i = 0; i < s; i++) {
				top = -1;
				int inner = get_inner_hit_id(doublets, i, layer);
				thisCell[0] = i;
				thisCell[1] = inner;
				thisCell[2] = get_outer_hit_id(doublets, i, layer);
				
				int optr = outerptr[layer-1][inner];
				while (optr != -1) {
					otherCell[0] = fetch_ml (isOuterHitOfCell[layer-1], optr);
					otherCell[1] = get_inner_hit_id(doublets, otherCell[0], layer-1);
					otherCell[2] = get_outer_hit_id(doublets, otherCell[0], layer-1);

					if (check_alignment_and_tag(doublets, thisCell, layer, otherCell,
								fargs[0], fargs[1], fargs[2],
								fargs[3], fargs[4], fargs[5])) {
						val = otherCell[0];
						top = push_back_ml (ml, top, val);
						if (top < 0) {
							printf ("Error: out of space\n");
							results.m_size = -1;
							flag = 1;
							break;
						}
					}
					optr = next_ml (isOuterHitOfCell[layer-1], optr);
				}
				liptr[layer][i] = top;
			}

			if (flag == 1)
				goto res;
		}
	

		unsigned int lastLayerPairIndex = numberOfLayers - 2;
	
		resetrv(&results);	

		OMPSimpleVector stack;

		/*get the quadruplets*/
		s = doublets[lastLayerPairIndex].size;
		#pragma omp parallel for num_threads(TC) private(stack)
		for (i = 0; i < s; i++) {
			resetsv(&stack);
			push_backsv(&stack, i);
			find_ntuplets(doublets, lastLayerPairIndex, i, ml, liptr, &results, &stack, 4);
		}

		res:
		if (results.m_size == -1)
			fail++;

		/*starts sending results and getting new parameters*/
		mppa_aiocb_set_trigger(d_portal_aiocb, 1);
		status |= mppa_aio_read(d_portal_aiocb);
		status |= mppa_pwrite(c_portal_fd, &results, sizeof(OMPResultVector), 0);
		mppa_aio_write(c_portal_aiocb);	
	}

	printf ("Failed to compute: %d\n", fail);
	

	mppa_exit(0);
	return 0;
}

