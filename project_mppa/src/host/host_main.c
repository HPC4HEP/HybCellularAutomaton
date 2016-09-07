/*! \file host_main.c this file has the functions used by the host process. 
*This process retrieves problems from files and sends them to the MPPA to be computed.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <time.h>

#include <omp.h>

#include "settings.h"

#include "OMPHitsAndDoublets.h"
#include "OMPResultVector.h"

#include <mppaipc.h>
#include <pcie.h>

//! \brief A function for loading a track seeding doublets structure from a file
/*! This function retrieves a track seeding doublets structure from a file. It also retrieves some floats that are used for compatibility 
* testing which are also included in the file.File name format is specified in settings.h (e.g. input_[var]_[var].txt)
\param x Identifier of the event (in)
\param y Identifier of the barrel (in)
\param fptr An array of 6 float used in the compatibility check afterwards (out)
\param iptr An array of sizes for each layer of the doublets structure that is sent to the compute cluster processing (out)
\param doublets The data structure for storing the doublets and hits of the problem (out)
*/
void dataLoad (int x, int y, float* fptr, int* iptr, OMPLayerDoublets* doublets) {
	char filename[100];

	float data;
	/*open input file*/
	sprintf (filename, inpath, x, y);
	printf ("File: %s\n", filename);
	FILE* fin = fopen(filename, "r");
	if (fin == NULL)
		printf ("Error!!!\n");

	int i;
	/*get float parameters*/
	for (i = 0; i < 6; i++) {
		fscanf (fin, "%f", &data);
		fptr[i] = data;		
	}
	/*retrieve the data structure*/
	int j,k;
	for (i = 0; i < theNumberOfLayers-1; i++) {
		fscanf (fin, "%d", &(doublets[i].size));
		iptr[3*i] = doublets[i].size;

		doublets[i].indices = (int*) malloc(2*doublets[i].size*sizeof(int));
		for (j = 0; j < 2*doublets[i].size; j++)
			fscanf (fin, "%d", &(doublets[i].indices[j]));
		
		for (k = 0; k < 2; k++) {
			fscanf (fin, "%d", &(doublets[i].layers[k].size));
			iptr[3*i+k+1] = doublets[i].layers[k].size;
			doublets[i].layers[k].p = (float*) malloc(3*doublets[i].layers[k].size*sizeof(float));
			
			for (j = 0; j < doublets[i].layers[k].size; j++) {
				fscanf (fin, "%f", &(doublets[i].layers[k].p[3*j]));
				fscanf (fin, "%f", &(doublets[i].layers[k].p[3*j+1]));
				fscanf (fin, "%f", &(doublets[i].layers[k].p[3*j+2]));
			}
		}
	}

	fclose (fin);

}

//! \brief The basic functionality of the host process is here. This function sents problems to clusters and gets results.
/*! This function gets file descriptors for communicating with I/O clusters and the problems and commits them to be solved. It then retrieves the results and prints the size. File descriptors are 2-element arrays, one for each I/O cluster.
\param fd_results The file descriptors for the buffers in which results are received
\param fd_l0 The file descriptor for the buffers from which layer 0 is sent
\param fd_l1 The file descriptor for the buffers from which layer 1 is sent
\param fd_l2 The file descriptor for the buffers from which layer 2 is sent
\param fd_args The file descriptor for the buffers from which the float parameters and sizes related to the hit and doublet structure is sent
\param fptr The arrays with the float parameters for all problems
\param iptr The arrays with the sizes related to the hit and doublet structures for all problems
\param doubletsList The hit and doublet structure for all problems
*/
void host_compute (int* fd_results, int* fd_l0, int* fd_l1, int* fd_l2, int* fd_args, float** fptr, int** iptr, OMPLayerDoublets** doubletsList) {
	printf ("host starts\n");
	int i,j,idx;
	int num[2][3];
	char* buf[2][3];

	/*
	This version is just a test for the given problem benchmark in mppa-256

	In a version that could be integrated to the full system, an idea would be to modify this in a client-server architecture as follows:
	-make mppa interface operating in a different thread and on a structure shared with client
	-each I/O cluster can have a queue for pending/finished tasks
		-possibly multithreaded I/O cluster with consumer-producer architecture so that it can enqueue more jobs
	-a function of host structure enqueues tasks (blocking or reject if full?) and gets task id
	-the independent mppa process submits jobs and retrieves results periodically by communicating with some thread of the I/O clusters
	-client retrieves results by task id
	*/

	int flag[2];
	flag[0] = 0;
	flag[1] = 0;
	/*memory allocated now so that we can have asynchronous communication*/
	buf[0][0] = malloc(1000000);
	buf[0][1] = malloc(1000000);
	buf[0][2] = malloc(1000000);
	buf[1][0] = malloc(1000000);
	buf[1][1] = malloc(1000000);
	buf[1][2] = malloc(1000000);

	OMPResultVector result[2];
	/*result async structures*/
	mppa_aiocb_t r_buffer_aiocb[2] = {	MPPA_AIOCB_INITIALIZER(fd_results[0], &result[0], sizeof(OMPResultVector)),
						MPPA_AIOCB_INITIALIZER(fd_results[1], &result[1], sizeof(OMPResultVector))};
	/*write async structures for parameter and layer sending*/
	mppa_aiocb_t wa_portal_aiocb[2];
	mppa_aiocb_t w1_portal_aiocb[2];
	mppa_aiocb_t w2_portal_aiocb[2];
	mppa_aiocb_t w3_portal_aiocb[2];

	mppa_aio_read(&r_buffer_aiocb[0]);
	mppa_aio_read(&r_buffer_aiocb[1]);

	for (idx = 0; idx < 700; idx++)
		{
			i = idx/7;
			j = idx%7;
			int n;
			/*switch between the clusters*/
			int u = (7*i+j)%2;
			OMPLayerDoublets* doublets = doubletsList[7*i+j];
			float* fargs = fptr[7*i+j];
			int* iargs = iptr[7*i+j];
			/*wait for previous writes to this I/O cluster to finish*/
			if (flag[u]) {
				mppa_aio_wait(&w1_portal_aiocb[u]);
				mppa_aio_wait(&w2_portal_aiocb[u]);
				mppa_aio_wait(&w3_portal_aiocb[u]);
			}
			/*re-using a cluster, so have to get result now*/
			if (7*i+j >= 16) {
				mppa_aio_wait(&r_buffer_aiocb[u]); 

				n = sizerv(&result[u]);
				printf ("%d:Size of results is %d\n", 7*i+j-16, n);

				mppa_aio_read(&r_buffer_aiocb[u]);

				char filename[100];

				int l = 7*i+j-16;

				int s1 = l/7;
				int s2 = l%7;

				sprintf (filename, outpath, s1, s2);
			}

			int num_a = 6*sizeof(float)+3*(theNumberOfLayers-1)*sizeof(int);

			char args[num_a];

			memcpy(args, fargs, 6*sizeof(float));
			memcpy(args+6*sizeof(float), iargs, 3*(theNumberOfLayers-1)*sizeof(int));

			num[u][0] = 2*sizeof(int)*doublets[0].size+3*sizeof(float)*(doublets[0].layers[0].size+doublets[0].layers[1].size);
			num[u][1] = 2*sizeof(int)*doublets[1].size+3*sizeof(float)*doublets[1].layers[1].size;
			num[u][2] = 2*sizeof(int)*doublets[2].size+3*sizeof(float)*doublets[2].layers[1].size;
			/*init asynchronous writes so that the they are ready*/
			mppa_aiocb_ctor(&w1_portal_aiocb[u], fd_l0[u], buf[u][0], num[u][0]);
			mppa_aiocb_set_pwrite(&w1_portal_aiocb[u], buf[u][0], num[u][0], 0);
			mppa_aiocb_ctor(&w2_portal_aiocb[u], fd_l1[u], buf[u][1], num[u][1]);
			mppa_aiocb_set_pwrite(&w2_portal_aiocb[u], buf[u][1], num[u][1], 0);
			mppa_aiocb_ctor(&w3_portal_aiocb[u], fd_l2[u], buf[u][2], num[u][2]);
			mppa_aiocb_set_pwrite(&w3_portal_aiocb[u], buf[u][2], num[u][2], 0);
			/*send parameters, then the I/O will start the reads for the asynchronous*/
			if (mppa_pwrite (fd_args[u], args, num_a, 0) != num_a) {
				printf ("failed to write args\n");
				exit(1);
			}

			int offset;
			int size;
			
			offset = 0;
			size = 2*sizeof(int)*doublets[0].size;
			memcpy (buf[u][0], doublets[0].indices, size);
			offset += size;
			size = 3*sizeof(float)*doublets[0].layers[0].size;
			memcpy (buf[u][0]+offset, doublets[0].layers[0].p, size);
			offset += size;
			size = 3*sizeof(float)*doublets[0].layers[1].size;
			memcpy (buf[u][0]+offset, doublets[0].layers[1].p, size);
			/*write layer 0*/
			mppa_aio_write(&w1_portal_aiocb[u]);

			offset = 0;
			size = 2*sizeof(int)*doublets[1].size;
			memcpy (buf[u][1], doublets[1].indices, size);
			offset += size;
			size = 3*sizeof(float)*doublets[1].layers[1].size;
			memcpy (buf[u][1]+offset, doublets[1].layers[1].p, size);
			/*write layer 1*/
			mppa_aio_write(&w2_portal_aiocb[u]);

			offset = 0;
			size = 2*sizeof(int)*doublets[2].size;
			memcpy (buf[u][2], doublets[2].indices, size);
			offset += size;
			size = 3*sizeof(float)*doublets[2].layers[1].size;
			memcpy (buf[u][2]+offset, doublets[2].layers[1].p, size);
			/*write layer 2*/
			mppa_aio_write(&w3_portal_aiocb[u]);

			flag[u] = 1;		
		}
	/*wait for the last writes to end*/
	mppa_aio_wait(&w1_portal_aiocb[0]);
	mppa_aio_wait(&w2_portal_aiocb[0]);
	mppa_aio_wait(&w3_portal_aiocb[0]);

	mppa_aio_wait(&w1_portal_aiocb[1]);
	mppa_aio_wait(&w2_portal_aiocb[1]);
	mppa_aio_wait(&w3_portal_aiocb[1]);
	
	/*retrieve last results*/
	for (i = 0; i < 16; i++) {
		int l = 700+i-16;

		int u = l%2;

		int s1 = l/7;
		int s2 = l%7;

		
		mppa_aio_wait(&r_buffer_aiocb[u]); 

		int n = sizerv(&result[u]);
		printf ("%d:Size of results is %d\n", l, n);

		char filename[100];

		sprintf (filename, outpath, s1, s2);

		if (i < 14)
			mppa_aio_read(&r_buffer_aiocb[u]); 

	}
}



int main(int argc, char *argv[]) {
	int load_ret[2];
	int mppa_status;
	int mppa_pid[2];

	char fargs[101];
	char iargs[101];
	
	memset(fargs, 0, 101);
	memset(iargs, 0, 101);

	const char* multibin_name = argv[1];
	const char* exec_name = argv[2];
	const char* cluster_exec_name = argv[3];

	int i;

	int j, k;

	float* fptr[700];
	int* iptr[700];
	OMPLayerDoublets* doublets[700];

	/*read problems from files*/
	for (i = 0; i < 100; i++)
		for (j = 0; j < 7; j++) {
			doublets[i*7+j] = malloc((theNumberOfLayers-1)*sizeof(OMPLayerDoublets));
			fptr[i*7+j] = malloc(6*sizeof(float));
			iptr[i*7+j] = malloc(3*(theNumberOfLayers-1)*sizeof(float));
			dataLoad (i, j, fptr[i*7+j], iptr[i*7+j], doublets[7*i+j]);
		}

	const char* io_args[3];
	io_args[0] = exec_name;
	io_args[1] = cluster_exec_name;
	io_args[2] = "0";
	io_args[3] = NULL;

	printf ("mppa host starts\n");
	
	struct stat buf;
	/*execute i/o processes*/
	if (stat(multibin_name, &buf)) {
		printf ("Error multi\n");
		exit(1);
	}

	if (stat(exec_name, &buf)) {
                printf ("Error exec\n");
                exit(1);
        }

	if (stat(cluster_exec_name, &buf)) {
                printf ("Error cluster\n");
                exit(1);
        }

	time_t start = time(NULL);

	if ((load_ret[0] = mppa_load(0,0,0,multibin_name)) < 0) {
		printf ("mppa load failed\n");
		exit(1);
	}


	if ((mppa_pid[0] = mppa_spawn(load_ret[0], NULL, exec_name, io_args, 0)) < 0) {
		printf ("mppa spawn failed\n");
		exit(1);
	}

	if ((load_ret[1] = mppa_load(0,0,1,multibin_name)) < 0) {
		printf ("mppa load failed\n");
		exit(1);
	}

	io_args[2] = "1";

	if ((mppa_pid[1] = mppa_spawn(load_ret[1], NULL, exec_name, io_args, 0)) < 0) {
		printf ("mppa spawn failed\n");
		exit(1);
	}

	/*open fds*/
	int fd_results[2];
	int fd_l0[2];
	int fd_l1[2];
	int fd_l2[2];
	int fd_args[2];

	fd_args[0] = mppa_open(args_name, O_WRONLY);
	if (fd_args[0] < 0) {
		printf ("failed to open args buffer\n");
		exit(1);
	}

	fd_l0[0] = mppa_open(l0_name, O_WRONLY);
	if (fd_l0[0] < 0) {
		printf ("failed to open l0 buffer\n");
		exit(1);
	}

	fd_l1[0] = mppa_open(l1_name, O_WRONLY);
        if (fd_l1[0] < 0) {
                printf ("failed	to open	l1 buffer\n");
                exit(1);
        }

	fd_l2[0] = mppa_open(l2_name, O_WRONLY);
        if (fd_l2[0] < 0) {
                printf ("failed	to open	l2 buffer\n");
                exit(1);
        }

	fd_results[0] = mppa_open(results_name, O_RDONLY);
        if (fd_results[0] < 0) {
                printf ("failed	to open	results	buffer\n");
                exit(1);
        }
 
	fd_args[1] = mppa_open(args_name2, O_WRONLY);
	if (fd_args[1] < 0) {
		printf ("failed to open args buffer\n");
		exit(1);
	}

	fd_l0[1] = mppa_open(l0_name2, O_WRONLY);
	if (fd_l0[1] < 0) {
		printf ("failed to open l0 buffer\n");
		exit(1);
	}

	fd_l1[1] = mppa_open(l1_name2, O_WRONLY);
        if (fd_l1[1] < 0) {
                printf ("failed	to open	l1 buffer\n");
                exit(1);
        }

	fd_l2[1] = mppa_open(l2_name2, O_WRONLY);
        if (fd_l2[1] < 0) {
                printf ("failed	to open	l2 buffer\n");
                exit(1);
        }

	fd_results[1] = mppa_open(results_name2, O_RDONLY);
        if (fd_results[1] < 0) {
                printf ("failed	to open	results	buffer\n");
                exit(1);
        }

		
	/*main work here*/	
	host_compute (fd_results, fd_l0, fd_l1, fd_l2, fd_args, fptr, iptr, doublets);

	time_t diff = time(NULL)-start;
	printf ("Execution time: %f sec\n", ((double) diff));

	/*unload i/o processes*/

	if ((mppa_waitpid(mppa_pid[0], &mppa_status, 0)) < 0) {
		printf ("mppa waitpid failed\n");
                exit(1);
	}

	if ((mppa_waitpid(mppa_pid[1], &mppa_status, 0)) < 0) {
		printf ("mppa waitpid failed\n");
                exit(1);
	}

	if ((mppa_unload(load_ret[0])) < 0) {
                printf ("mppa waitpid failed\n");
                exit(1);
        }

	if ((mppa_unload(load_ret[1])) < 0) {
                printf ("mppa waitpid failed\n");
                exit(1);
        }
	
	return 0;
}
