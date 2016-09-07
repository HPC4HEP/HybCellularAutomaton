#ifndef __SHARED_DEFS_H__
#define __SHARED_DEFS_H__

#define IO_EXECUTABLE "project_io"
#define CLUSTER_EXECUTABLE "project_cluster"

#define CLUSTER_COUNT 	16
#define CLUSTER_RANGE 	"[0..15]"

#define a_buffer_name	"/mppa/buffer/board0#mppa0#pcie0#0/host#0"
#define b_buffer_name	"/mppa/buffer/board0#mppa0#pcie0#1/host#1"
#define c_buffer_name	"/mppa/buffer/host#2/board0#mppa0#pcie0#2"

/**
 * Vector length shared between host and IO
 * This could easily be dynamic
 */
#define VECTOR_LENGTH	100000

#ifdef __K1__
#  ifdef __k1io__
#    define ARCH	"io"
#  else
#    define ARCH	"cluster"
#  endif
#else
#  define ARCH	"host"
#endif

#define test_printf(__fmt, ...) printf("["ARCH" test] "__fmt, ##__VA_ARGS__)

#endif /* __SHARED_DEFS_H__ */
