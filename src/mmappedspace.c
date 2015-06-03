/*
 * A simple dlmalloc wrapper
 *
 * Kazutomo Yoshii <ky@anl.gov>
 */

#include <unistd.h>
#include <sys/types.h>

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <error.h>
#include <string.h>

#include "mmappedspace.h"

/* dlmalloc definitions.  */
typedef void *mspace;
extern mspace create_mspace_with_base(void* base, size_t capacity, int locked);
extern void *mspace_malloc(mspace msp, size_t bytes);
extern void  mspace_free(mspace msp, void* mem);

/*---*/

void *mms_malloc(struct mms_struct *mms, size_t size)
{
	if (!mms) {
		return NULL;
	}
	return mspace_malloc(mms->mspace_ptr, size);
}

void mms_free(struct mms_struct *mms, void* ptr)
{
	if (!mms) {
		return;
	}
	return mspace_free(mms->mspace_ptr, ptr);
}


void *mms_get_ptr(struct mms_struct *mms)
{
	if (!mms) {
		return NULL;
	}
	return mms->ptr;
}

struct mms_struct *mms_init(const char *fn, size_t totalsize)
{
	struct mms_struct *mms;

	mms = malloc(sizeof(struct mms_struct));
	if (!mms) {
		return NULL;
	}

	mms->size = totalsize;

	snprintf(mms->fn, MMS_FN_MAXLEN, "%s", fn);
#if 0	
	if (access(mms->fn, F_OK) == 0) {
		printf("Remove %s first\n", mms->fn);
		goto mms_free;
	}
#endif
	mms->fd = open(mms->fn, O_CREAT|O_RDWR, 0755);
	if (mms->fd < 0) {
		perror("open");
		goto mms_free;
	}
#if 0
	if (ftruncate(mms->fd, mms->size) ) {
		perror("ftruncate");
		goto mms_unlink;
	}
#endif

	mms->ptr = mmap(0, mms->size, PROT_READ|PROT_WRITE, 
			MAP_SHARED, mms->fd, 0);
	if (mms->ptr == MAP_FAILED) {
		perror("mmap");
		close(mms->fd);
		goto mms_unlink;
	}

	mms->mspace_ptr = create_mspace_with_base(mms->ptr, mms->size, 0);

	if (!mms->mspace_ptr) {
		perror("create_mspace_with_base");
		goto mms_unmap;
	}

	return mms;

 mms_unmap:
	munmap(mms->ptr, mms->size);
 mms_unlink:
//	unlink(mms->fn);
 mms_free:
	free(mms);

	return NULL;
}

void  mms_fini(struct mms_struct *mms)
{
	if (!mms) {
		return ;
	}

	if (mms->ptr != MAP_FAILED) {
		munmap(mms->ptr, mms->size);
		//unlink(mms->fn);
	}
}

#ifdef __MMAPPEDSPACE_TEST_MAIN__

int main(int argc, char *argv[])
{
	size_t totalsize;
	size_t testsize;
	char *fn;
	int  size_shift;
	struct mms_struct *mms;
	unsigned long *allocptrs[4];
	int i;

	if (argc < 3) {
		printf("%s fn size_shift\n",
			argv[0]);
		exit(1);
	}

	fn = argv[1];
	size_shift = atoi(argv[2]);

	if (size_shift<16) {
		printf("size_shift should be >= 16!\n");
		exit(1);
	}

	totalsize = 1ULL << size_shift;
	testsize =  totalsize >> 2;

	printf("totalsize = %lu\n", totalsize);
	printf("testsize  = %lu\n", testsize);

	mms = mms_init(fn, totalsize);
	if (!mms) {
		printf("mms_init() failed\n");
		exit(1);
	}

	for(i = 0; i < 4 ; i++) {
		allocptrs[i] = mms_malloc(mms, testsize);
		if (!allocptrs[i]) {
			printf("Failed to allocate %u bytes\n", 1<<i);
			exit(1);
		}
		printf("%2d : %p allocated\n", i, allocptrs[i]);
	}

	for(i = 0; i < 4 ; i++) {
		memset(allocptrs[i], 0, testsize);
		printf("%2d : %p memset\n", i, allocptrs[i]);
	}

	for(i = 0; i < 4 ; i++) {
		mms_free(mms, allocptrs[i]);
		printf("%2d : %p free'ed\n", i, allocptrs[i]);
	}

	mms_fini(mms);

	return 0;
}

#endif
