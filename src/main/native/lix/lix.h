#ifndef LINEAR_INDEX_LIX_H
#define LINEAR_INDEX_LIX_H

#include <stdint.h>
#include "htslib/tbx.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"


#ifdef __cplusplus
extern "C" {
#endif

#define LIX_VERSION "0.2.2"


typedef struct linear_offset {
    int64_t n, m;
    uint64_t *offset;
} linear_offset_t;


typedef struct lix {
    int32_t fmt, min_shift;
    int32_t n;
    linear_offset_t *linear_offsets;
    void *dict;
} lix_t;


/**
 * Create lix_t from tbx_t by copying linear index from tbx_t, 
 * min_shift is the minimum shift which specifies the interval size(1<<min_shift) in lix_t.
 * 
 * @param tbx htslib tbx_t
 * @param min_shift minimum shift, the interval size is 1 << min_shift
 * @return The lix_t index, or NULL if an error occurred
 */
lix_t *lix_create(const tbx_t *tbx, int32_t min_shift);


/**
 * Free lix memory.
 * @param lix The lix index.
 */
void lix_destroy(lix_t *lix);



/**
 * Save the lix index to file.
 * @param lix The lix index.
 * @param lix_fn Output file name.
 * @return 0 (success), < 0 (failure)
 */
int lix_save(const lix_t *lix, const char *lix_fn);


/**
 * Build the lix index from tbi index file.
 * This function just calls tbx_index_load and calls lix_create to create the index.
 * @param fn file name of tbi index
 * @param min_shift minimum shift, the interval size is 1 << min_shift
 * @param out_fn output file name
 * @return 0 (success), -1 (failure)
 */
int lix_build(const char *fn, int32_t min_shift, const char *out_fn);


/**
 * Load the lix index.
 * @param lix_fn file name of lix index
 * @return The lix_t index, or NULL if an error occurred
 */
lix_t *lix_load(const char *lix_fn);


/**
 * Get names from lix
 * @param lix The lix index
 * @param name Output names, the memory should be allocated by the caller
 * @return number of names (success), < 0 (failure)
 */
int lix_get_names(const lix_t *lix, const char **names);


/**
 * Get id(rank) of the name from the index
 * @param lix The lix index
 * @param name name
 * @return The id, or -1 if name is not in the index
 */
int32_t lix_name2id(const lix_t *lix, const char *name);


/**
 * Check if two lix index is equal
 * @param lix_a The first lix index
 * @param lix_b The second lix index
 * @return 0 (equal), 1 (not equal), -1 (failure)
 */
int lix_equal(const lix_t *lix_a, const lix_t *lix_b);


typedef struct lix_iter {
    int32_t has_next;
    char *contig;
    int64_t start, end;
    uint64_t offset;
} lix_iter_t;

/**
 * Initialize lix_iter_t 
 */
lix_iter_t *lix_iter_init(void);

/**
 * Free lix_iter_t memory
 * @param iter iterator
 */
void lix_iter_destroy(lix_iter_t *iter);

/**
 * Seek to the offset of the given region
 * @param idx lix index
 * @param contig query contig name
 * @param start 0-based query start
 * @param end 0-based query end
 * @return iterator, the memory should be freed by lix_iter_desdroy
 */
lix_iter_t *lix_query(const lix_t *idx, const char *contig, int64_t start, int64_t end);

/**
 * Read next record from the iter/fp
 * @param iter iterator
 * @param fp vcf file pointer
 * @param hdr vcf file header
 * @param rec next record
 * @return 0 on success, -1 if there is no more records, < -1 on error
 */
int lix_iter_next(lix_iter_t *iter, htsFile *fp, const bcf_hdr_t *hdr, bcf1_t *rec);


#ifdef __cplusplus
}
#endif


#endif  // LINEAR_INDEX_LIX_H

