/*! \file io_main.c this file has the functions used by the i/o cluster process. 
*This process retrieves problems from the host and forwards them to clusters. Afterwards, it retrieves the results.
*/

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "shared_defs.h"
#include "settings.h"


#include "OMPResultVector.h"
#include "OMPHitsAndDoublets.h"
#include "OMPCACell.h"

#include <mppaipc.h>
#include <mppa/osconfig.h>
#include <pcie_queue.h>


/*buffers, one for each cluster so that we can support asynchronous communication and pipelined execution*/
char* vector_l0[8];
char* vector_l1[8];
char* vector_l2[8];
OMPResultVector* vector_c[8];
const int num_a = 6*sizeof(float)+3*3*sizeof(int);

//! \brief The basic functionality of the io process is here.
/*! This function gets file descriptors for communicating with compute clusters and the host, and schedules problem solving. In each iteration,
* its choose the cluster and transfers any solution it might have to the host, then
* it gets info on the size of received problems and receives the data. As soon as a portion of data is received (layer pair),
* start forwarding it to the compute clusters.
\param io_id The identifier of the io cluster so that we can point to the correct resources
\param fd_host An array of file descriptors for communicating with the host
\param fd_cluster An array of file descriptors for communicating with the clusters
*/
int io_compute (int io_id, int* fd_host, int** fd_cluster) {
	int i;
	int num_i, num_p;
	int length;
	int l1, l2;
	
	char vector_d[num_a];	
	int err;

	int num[3];
	int iargs[3*3];
	/*allocate buffers for asynchronous communication*/
	for (i = 0; i < 8; i++) {
		vector_l0[i] = malloc(1000000);
		vector_l1[i] = malloc(1000000);
		vector_l2[i] = malloc(1000000);
		vector_c[i] = malloc(sizeof(OMPResultVector));
	}
	/*set up asynchronous channels*/
	mppa_aiocb_t r_portal_aiocb[8] =
				{ MPPA_AIOCB_INITIALIZER(fd_cluster[4][io_id*8+0], vector_c[0], sizeof(OMPResultVector)),
				  MPPA_AIOCB_INITIALIZER(fd_cluster[4][io_id*8+1], vector_c[1], sizeof(OMPResultVector)),
				  MPPA_AIOCB_INITIALIZER(fd_cluster[4][io_id*8+2], vector_c[2], sizeof(OMPResultVector)),
				  MPPA_AIOCB_INITIALIZER(fd_cluster[4][io_id*8+3], vector_c[3], sizeof(OMPResultVector)),
				  MPPA_AIOCB_INITIALIZER(fd_cluster[4][io_id*8+4], vector_c[4], sizeof(OMPResultVector)),
				  MPPA_AIOCB_INITIALIZER(fd_cluster[4][io_id*8+5], vector_c[5], sizeof(OMPResultVector)),
				  MPPA_AIOCB_INITIALIZER(fd_cluster[4][io_id*8+6], vector_c[6], sizeof(OMPResultVector)),
				  MPPA_AIOCB_INITIALIZER(fd_cluster[4][io_id*8+7], vector_c[7], sizeof(OMPResultVector))};
	/*async communication structures for writing to clusters*/
	mppa_aiocb_t w1_portal_aiocb[8];
	mppa_aiocb_t w2_portal_aiocb[8];
	mppa_aiocb_t w3_portal_aiocb[8];


	for (l1 = 0; l1 < 350; l1+=8) {
		int step = (l1 + 8 > 350)? 350-l1: 8;

		for (l2 = 0; l2 < step; l2++) {
			int cluster = io_id*8+l2;
			/*fetch result from cluster that is about to be reused*/
			if (l1 >= 8) {
				mppa_aio_wait(&w1_portal_aiocb[cluster%8]);
				mppa_aio_wait(&w2_portal_aiocb[cluster%8]);
				mppa_aio_wait(&w3_portal_aiocb[cluster%8]);

				mppa_aio_wait(&r_portal_aiocb[l2]);
				length = sizeof(OMPResultVector);
				if(mppa_pwrite(fd_host[2], vector_c[l2], length, 0) != length) {
					test_printf("mppa_write fd_c failed\n");
					exit(1);
				}
			}
			/*get arguments*/
			if((err = mppa_read(fd_host[3], vector_d, num_a)) != num_a) {
				test_printf("mppa_read fd_d failed %d\n", err);
				exit(1);
			}

			memcpy(iargs, vector_d+6*sizeof(float), 3*3*sizeof(int));

			num[0] = 2*sizeof(int)*iargs[0]+3*sizeof(float)*(iargs[1]+iargs[2]);
			num[1] = 2*sizeof(int)*iargs[3]+3*sizeof(float)*iargs[5];
			num[2] = 2*sizeof(int)*iargs[6]+3*sizeof(float)*iargs[8];

			/*initialize receives*/
			mppa_aiocb_t l0_portal_aiocb[1] =
				{ MPPA_AIOCB_INITIALIZER(fd_host[0], vector_l0[cluster%8], num[0]) };
			mppa_aiocb_t l1_portal_aiocb[1] =
				{ MPPA_AIOCB_INITIALIZER(fd_host[1], vector_l1[cluster%8], num[1]) };
			mppa_aiocb_t l2_portal_aiocb[1] =
				{ MPPA_AIOCB_INITIALIZER(fd_host[4], vector_l2[cluster%8], num[2]) };
			
			mppa_aio_read(l0_portal_aiocb);
			mppa_aio_read(l1_portal_aiocb);
			mppa_aio_read(l2_portal_aiocb);


			mppa_aiocb_ctor(&w1_portal_aiocb[cluster%8], fd_cluster[2][cluster], vector_l0[cluster%8], num[0]);
			mppa_aiocb_set_pwrite(&w1_portal_aiocb[cluster%8], vector_l0[cluster%8], num[0], 0);
			mppa_aiocb_ctor(&w2_portal_aiocb[cluster%8], fd_cluster[3][cluster], vector_l1[cluster%8], num[1]);
			mppa_aiocb_set_pwrite(&w2_portal_aiocb[cluster%8], vector_l1[cluster%8], num[1], 0);
			mppa_aiocb_ctor(&w3_portal_aiocb[cluster%8], fd_cluster[5][cluster], vector_l2[cluster%8], num[2]);
			mppa_aiocb_set_pwrite(&w3_portal_aiocb[cluster%8], vector_l2[cluster%8], num[2], 0);
			
			long long dummy;

			mppa_aio_read(&r_portal_aiocb[l2]);
			/*sync ensures that communication will be initialized in both sides*/
			long long mask = -2;
			mppa_ioctl (fd_cluster[0][cluster], MPPA_RX_SET_MATCH, mask);
	
			if (mppa_read (fd_cluster[0][cluster], &dummy, sizeof(long long)) < 0) {
				printf ("Sync failed\n");
				return -1;
			}

			if (mppa_pwrite (fd_cluster[1][cluster], vector_d, num_a, 0) != num_a) {
				printf ("args failed\n");
				return -1;
			}

			mppa_ioctl (fd_cluster[0][cluster], MPPA_RX_SET_MATCH, mask);

			if (mppa_read (fd_cluster[0][cluster], &dummy, sizeof(long long)) < 0) {
				printf ("Sync failed\n");
				return -1;
			}
			/*forward data*/
			mppa_aio_wait(l0_portal_aiocb);
			mppa_aio_write(&w1_portal_aiocb[cluster%8]);

			mppa_aio_wait(l1_portal_aiocb);
			mppa_aio_write(&w2_portal_aiocb[cluster%8]);

			mppa_aio_wait(l2_portal_aiocb);
			mppa_aio_write(&w3_portal_aiocb[cluster%8]);
		}
	}

	for (i = 0; i < 8; i++) {
		int cluster = (i+350)%8+io_id*8;
		mppa_aio_wait(&r_portal_aiocb[cluster%8]);

		length = sizeof(OMPResultVector);
	
		if(mppa_pwrite(fd_host[2], vector_c[cluster%8], length, 0) != length) {
			test_printf("mppa_write fd_c failed\n");
			exit(1);
		}
	}

	return 0;
}


int main(int argc, char *argv[])
{
	int i, j;
	int iargs[3*3];

	int io_id = atoi(argv[2]);


	const char* root_sync = (io_id == 0)?"/mppa/sync/128:1":"/mppa/sync/192:1";
	const char* l0_portal[16] = {"/mppa/portal/0:2","/mppa/portal/1:2","/mppa/portal/2:2","/mppa/portal/3:2",
								"/mppa/portal/4:2","/mppa/portal/5:2","/mppa/portal/6:2","/mppa/portal/7:2",
								"/mppa/portal/8:2","/mppa/portal/9:2","/mppa/portal/10:2","/mppa/portal/11:2",
								"/mppa/portal/12:2","/mppa/portal/13:2","/mppa/portal/14:2","/mppa/portal/15:2"};
	const char* l1_portal[16] = {"/mppa/portal/0:5","/mppa/portal/1:5","/mppa/portal/2:5","/mppa/portal/3:5",
								"/mppa/portal/4:5","/mppa/portal/5:5","/mppa/portal/6:5","/mppa/portal/7:5",
								"/mppa/portal/8:5","/mppa/portal/9:5","/mppa/portal/10:5","/mppa/portal/11:5",
								"/mppa/portal/12:5","/mppa/portal/13:5","/mppa/portal/14:5","/mppa/portal/15:5"};
	const char* l2_portal[16] = {	"/mppa/portal/0:7","/mppa/portal/1:7","/mppa/portal/2:7","/mppa/portal/3:7",
					"/mppa/portal/4:7","/mppa/portal/5:7","/mppa/portal/6:7","/mppa/portal/7:7",
					"/mppa/portal/8:7","/mppa/portal/9:7","/mppa/portal/10:7","/mppa/portal/11:7",
					"/mppa/portal/12:7","/mppa/portal/13:7","/mppa/portal/14:7","/mppa/portal/15:7"};
	const char* a_portal[16] = {"/mppa/portal/0:8","/mppa/portal/1:8","/mppa/portal/2:8","/mppa/portal/3:8",
								"/mppa/portal/4:8","/mppa/portal/5:8","/mppa/portal/6:8","/mppa/portal/7:8",
								"/mppa/portal/8:8","/mppa/portal/9:8","/mppa/portal/10:8","/mppa/portal/11:8",
								"/mppa/portal/12:8","/mppa/portal/13:8","/mppa/portal/14:8","/mppa/portal/15:8"};
	const char* r_portal[16] = {"/mppa/portal/128:17","/mppa/portal/128:2","/mppa/portal/128:3","/mppa/portal/128:4",
					"/mppa/portal/128:5","/mppa/portal/128:6","/mppa/portal/128:7","/mppa/portal/128:8",
					"/mppa/portal/192:17","/mppa/portal/192:2","/mppa/portal/192:3","/mppa/portal/192:4",
					"/mppa/portal/192:5","/mppa/portal/192:6","/mppa/portal/192:7","/mppa/portal/192:8"};
	
	int pid[8];

	int k;
	/*spawn cluster processes*/
	for (k = 0; k < 8; k++) {
		i = io_id*8+k;
		const char* args[] = { "project_cluster", root_sync, a_portal[i], l0_portal[i], l1_portal[i], r_portal[i], l2_portal[i], NULL };
		if ((pid[k] = mppa_spawn (i, NULL, "project_cluster", args, 0)) < 0) {
			printf ("Spawn failed\n");
			return -1;
		}
	}

	int fd_a, fd_b, fd_c, fd_d, fd_e;
	size_t length;

	int fd_host[5];
	int** fd_cluster = malloc(6*sizeof(int*));
	for (i = 0; i < 6; i++)
		fd_cluster[i] = malloc(16*sizeof(int));

	const char *l0name, *l1name, *l2name, *rname, *aname;
	if (io_id == 0) {
		l0name = l0_name;
		l1name = l1_name;
		l2name = l2_name;
		rname = results_name;
		aname = args_name;
	} else {
		l0name = l0_name2;
		l1name = l1_name2;
		l2name = l2_name2;
		rname = results_name2;
		aname = args_name2;
	}
	/*open fds*/
	fd_a = mppa_open(l0name, O_RDONLY);
	if(fd_a < 0) {
		test_printf("open buffer l0 failed\n");
		exit(1);
	}
	fd_b = mppa_open(l1name, O_RDONLY);
	if(fd_b < 0) {
		test_printf("open buffer l1 failed\n");
		exit(1);
	}
	fd_c = mppa_open(rname, O_WRONLY);
	if(fd_c < 0) {
		test_printf("open buffer result failed\n");
		exit(1);
	}
	fd_d = mppa_open(aname, O_RDONLY);
	if(fd_d < 0) {
		test_printf("open buffer arg failed\n");
		exit(1);
	}
	fd_e = mppa_open(l2name, O_RDONLY);
	if(fd_e < 0) {
		test_printf("open buffer l2 failed\n");
		exit(1);
	}

	fd_host[0] = fd_a; 
	fd_host[1] = fd_b; 
	fd_host[2] = fd_c; 
	fd_host[3] = fd_d;
	fd_host[4] = fd_e;


	int root_sync_fd, a_portal_fd, l0_portal_fd,
		l1_portal_fd, l2_portal_fd, r_portal_fd;

	root_sync_fd = mppa_open(root_sync, O_RDONLY);

	for (k = 0; k < 8; k++) {
		i = io_id*8+k;
		a_portal_fd = mppa_open(a_portal[i], O_WRONLY);
		l0_portal_fd = mppa_open(l0_portal[i], O_WRONLY);
		l1_portal_fd = mppa_open(l1_portal[i], O_WRONLY);
		l2_portal_fd = mppa_open(l2_portal[i], O_WRONLY);
		r_portal_fd = mppa_open(r_portal[i], O_RDONLY);	

		fd_cluster[0][i] = root_sync_fd; 
		fd_cluster[1][i] = a_portal_fd; 
		fd_cluster[2][i] = l0_portal_fd; 
		fd_cluster[3][i] = l1_portal_fd; 
		fd_cluster[4][i] = r_portal_fd;
		fd_cluster[5][i] = l2_portal_fd;
	}

	/* perform the work */

	io_compute (io_id, fd_host, fd_cluster);
	
	return 0;
}
