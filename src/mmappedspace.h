/*
 * A simple dlmalloc wrapper
 *
 * Kazutomo Yoshii <ky@anl.gov>
 */
#ifndef __MMAPPEDSPACE_H_DEFINED__
#define __MMAPPEDSPACE_H_DEFINED__

#define  MMS_FN_MAXLEN  (256)
struct mms_struct {
	int    fd;
	size_t size;  /* size for mmap() */
	void   *ptr;
	void   *mspace_ptr;
	char   fn[MMS_FN_MAXLEN];
};


extern struct mms_struct *mms_init(const char *fn, size_t totalsize);
extern void  mms_fini(struct mms_struct *mms);

extern void *mms_malloc(struct mms_struct *mms, size_t size);
extern void mms_free(struct mms_struct *mms, void* ptr);




#endif
