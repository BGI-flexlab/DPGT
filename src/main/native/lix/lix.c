#include <unistd.h>
#include <string.h>
#include "lix.h"
#include "htslib/bgzf.h"
#include "htslib/hts_log.h"

// declare string to int64_t khash table
#include "htslib/khash.h"
// KHASH_DECLARE(s2i, kh_cstr_t, int64_t)
KHASH_MAP_INIT_STR(s2i, int64_t)

#define TBI_IDX_SUFFIX ".tbi"
#define CSI_IDX_SUFFIX ".csi"
#define LIX_IDX_SUFFIX ".lix"
#define LIX_FILE_MAGIC "LIX\1"


// copy from hst.c to export hts_idx_t struct
typedef struct {
    int32_t m, n;
    uint64_t loff;
    hts_pair64_t *list;
} bins_t;

KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
    hts_pos_t n, m;
    uint64_t *offset;
} lidx_t;

struct hts_idx_t {
    int fmt, min_shift, n_lvls, n_bins;
    uint32_t l_meta;
    int32_t n, m;
    uint64_t n_no_coor;
    bidx_t **bidx;
    lidx_t *lidx;
    uint8_t *meta; // MUST have a terminating NUL on the end
    int tbi_n, last_tbi_tid;
    struct {
        uint32_t last_bin, save_bin;
        hts_pos_t last_coor;
        int last_tid, save_tid, finished;
        uint64_t last_off, save_off;
        uint64_t off_beg, off_end;
        uint64_t n_mapped, n_unmapped;
    } z; // keep internal states
};
// end copy from hst.c


lix_t *lix_create(const tbx_t *tbx, int32_t min_shift)
{
    if (tbx == NULL) return NULL;
    lix_t *lix = (lix_t *)calloc(1, sizeof(lix_t));
    if (lix == NULL) {
        hts_log_error("Out of memory!");
        return NULL;
    }
    lix->fmt = tbx->idx->fmt;
    if (min_shift < 14) {
        // interval size is at least equal to 16K(1<<14), which is tbi's linear interval size
        hts_log_warning("min_shift is less than 14, set it to 14.");
        min_shift = 14;
    }
    lix->min_shift = min_shift;
    lix->n = tbx->idx->n;
    
    // copy linear index
    lix->linear_offsets = (linear_offset_t *)calloc(
        lix->n, sizeof(linear_offset_t));
    if (lix->linear_offsets == NULL) {
        hts_log_error("Out of memory!");
        free(lix);
        return NULL;
    }
    int fold = 1 << (min_shift - 14);
    int i, j;
    for (i = 0; i < lix->n; ++i) {
        lix->linear_offsets[i].n = (tbx->idx->lidx[i].n + fold - 1) / fold;
        lix->linear_offsets[i].m = lix->linear_offsets[i].n;
        lix->linear_offsets[i].offset = (uint64_t *)calloc(
            lix->linear_offsets[i].n, sizeof(uint64_t));
        if (lix->linear_offsets[i].offset == NULL) {
            hts_log_error("Out of memory!");
            for (j = 0; j < i; ++j) {
                free(lix->linear_offsets[j].offset);
            }
            free(lix->linear_offsets);
            free(lix);
            return NULL;
        }
        int k = 0;
        for (j = 0; j < tbx->idx->lidx[i].n; j += fold) {
            lix->linear_offsets[i].offset[k++] = tbx->idx->lidx[i].offset[j];
        }
    }

    lix->dict = kh_init(s2i);
    if (!lix->dict) {
        hts_log_error("Out of memory!");
        lix_destroy(lix);
        return NULL;
    }

    khint_t k, k1;
    int absent;
    khash_t(s2i) *d = (khash_t(s2i)*)tbx->dict;
    khash_t(s2i) *d1 = (khash_t(s2i)*)lix->dict;
    for (k = kh_begin(d); k != kh_end(d); ++k) {
        if (kh_exist(d, k)) {
            const char *name = kh_key(d, k);
            k1 = kh_put(s2i, d1, name, &absent);
            if (absent < 0) {
                hts_log_error("Out of memory!");
                lix_destroy(lix);
                return NULL;
            } else if (absent) {
                char *name_dup = strdup(name);
                if (name_dup) {
                    kh_key(d1, k1) = name_dup;
                    kh_val(d1, k1) = kh_val(d, k);
                } else {
                    kh_del(s2i, d1, k1);
                    hts_log_error("Out of memory!");
                    lix_destroy(lix);
                    return NULL; 
                }
            }
        }
    }

    return lix;
}


void lix_destroy(lix_t *lix)
{
    if (!lix) return;
    
    if (lix->dict) {
        khint_t k;
        khash_t(s2i) *d = (khash_t(s2i)*)lix->dict;
        // free keys
        for (k = kh_begin(d); k != kh_end(d); ++k) {
            if (kh_exist(d, k)) free((char*)kh_key(d, k));
        }
        kh_destroy(s2i, (khash_t(s2i) *)lix->dict);
    }
    lix->dict = 0;
    if (lix->linear_offsets) {
        int i;
        for (i = 0; i < lix->n; ++i) {
            free(lix->linear_offsets[i].offset);
        }
        free(lix->linear_offsets);
    }
    free(lix);
}


int lix_get_names(const lix_t *lix, const char **names)
{
    khint_t k;
    int n = 0;
    khash_t(s2i) *d = (khash_t(s2i)*)lix->dict;
    for (k = kh_begin(d); k != kh_end(d); ++k) {
        if (kh_exist(d, k)) {
            ++n;
            if (n > lix->n) return -1;
            names[kh_val(d, k)] = kh_key(d, k);
        }
    }
    if (n != lix->n) return -1;
    return n;
}


int32_t lix_name2id(const lix_t *lix, const char *name)
{
    khash_t(s2i) *d = (khash_t(s2i)*)lix->dict;
    if (!d) return -1;
    khint_t k = kh_get(s2i, d, name);
    if (k == kh_end(d)) {
        return -1;
    } else {
        return kh_val(d, k);
    }
}


int lix_save(const lix_t *lix, const char *lix_fn)
{
    #define check(ret) if ((ret) < 0) return -1

    BGZF *lix_fp = bgzf_open(lix_fn, "w");
    if (!lix_fp) {
        hts_log_error("Can not open %s for writing.", lix_fn);
        return -1;
    }

    check(bgzf_write(lix_fp, LIX_FILE_MAGIC, 4));
    check(bgzf_write(lix_fp, &lix->fmt, sizeof(lix->fmt)));
    check(bgzf_write(lix_fp, &lix->min_shift, sizeof(lix->min_shift)));
    check(bgzf_write(lix_fp, &lix->n, sizeof(lix->n)));
    
    const char **names = (const char **)malloc(lix->n*sizeof(const char *));
    int ret = lix_get_names(lix, names);
    if (ret < 0) {
        hts_log_error("Failed to get names from lix.");
        return -1;
    }
 
    int32_t length_of_names = 0;
    int i;
    for (i = 0; i < lix->n; ++i) {
        length_of_names += strlen(names[i]) + 1;
    }
    check(bgzf_write(lix_fp, &length_of_names, sizeof(length_of_names)));
    for (i = 0; i < lix->n; ++i) {
        check(bgzf_write(lix_fp, names[i], strlen(names[i])));
        check(bgzf_write(lix_fp, "\0", 1));
    }

    for (i = 0; i < lix->n; ++i) {
        check(bgzf_write(lix_fp, &lix->linear_offsets[i].n, sizeof(int64_t)));
        check(bgzf_write(lix_fp, lix->linear_offsets[i].offset,
            lix->linear_offsets[i].n*sizeof(uint64_t)));
    }

    free(names);
    ret = bgzf_close(lix_fp);
    if (ret < 0) return -1;
    return 0;

    #undef check
}


static int ends_with_tbi(const char *str) {
    const char *suffix = ".tbi";
    size_t str_len = strlen(str);
    size_t suffix_len = strlen(suffix);

    if (str_len >= suffix_len &&
        strcmp(str + str_len - suffix_len, suffix) == 0)
    {
        return 1;
    } else {
        return 0;
    }
}


int lix_build(const char *fn, int32_t min_shift, const char *out_fn)
{
    // only support .tbi index file
    if (!ends_with_tbi(fn)) {
        hts_log_error("Only .tbi index file is supported.");
        return -1;
    }

    tbx_t *tbx = tbx_index_load2(fn, fn);
    if (tbx == NULL) {
        return -1;
    }

    lix_t *lix = lix_create(tbx, min_shift);
    if (lix == NULL) {
        return -1;
    }

    int ret = lix_save(lix, out_fn);
    return ret;
}


static inline int check_bgzf_read(int ret, int n) {
    if (ret < 0 || ret != n) return -1;
    return 0;
}


lix_t *lix_load(const char *lix_fn)
{
    char *names = NULL;

    if (access(lix_fn, F_OK))
    {
        hts_log_error("lix file %s is not exists.", lix_fn);
        return NULL;
    }

    if (access(lix_fn, R_OK)) {
        hts_log_error("Permission denied, can not read %s %s%s", lix_fn,
            errno ? "error message: " : "", strerror(errno));
        return NULL;
    }

    BGZF *lix_fp = bgzf_open(lix_fn, "r");
    if (!lix_fp) {
        hts_log_error("Can not read %s %s%s", lix_fn,
            errno ? "error message: " : "", strerror(errno));
        return NULL;
    }

    char magic[4];
    int r = bgzf_read(lix_fp, magic, 4);
    if (r != 4) {
        hts_log_error("%s lix file format error! Can not read magic from it.",
            lix_fn);
        bgzf_close(lix_fp);
        return NULL;
    } else if (memcmp(magic, LIX_FILE_MAGIC, 4)) {
        char magic_str[5];
        strcat(magic_str, LIX_FILE_MAGIC);
        magic_str[4] = '\0';
        hts_log_error("%s lix file format error! File magic not match %s",
            lix_fn, magic_str);
        bgzf_close(lix_fp);
        return NULL;
    }

    lix_t *lix = (lix_t *)calloc(1, sizeof(lix_t));
    if (!lix) hts_log_error("Out of memory!");

    int ret = check_bgzf_read(
        bgzf_read(lix_fp, &lix->fmt, sizeof(lix->fmt)),
        sizeof(lix->fmt));
    if (ret < 0) goto fail;
    ret = check_bgzf_read(
        bgzf_read(lix_fp, &lix->min_shift, sizeof(lix->min_shift)),
        sizeof(lix->min_shift));
    if (ret < 0) goto fail;
    ret = check_bgzf_read(
        bgzf_read(lix_fp, &lix->n, sizeof(lix->n)),
        sizeof(lix->n));
    if (ret < 0) goto fail;

    int32_t length_of_names = 0;
    ret = check_bgzf_read(
        bgzf_read(lix_fp, &length_of_names, sizeof(length_of_names)),
        sizeof(length_of_names));
    if (ret < 0) goto fail;
    
    names = (char *)malloc(length_of_names);
    if (!names) {
        hts_log_error("Out of memory!");
        goto fail;
    }

    ret = check_bgzf_read(
        bgzf_read(lix_fp, names, length_of_names),
        length_of_names);
    if (ret < 0) goto fail;

    lix->dict = kh_init(s2i);
    if (!lix->dict) {
        hts_log_error("Out of memory!");
        goto fail;
    }

    khint_t k;
    int absent;
    khash_t(s2i) *d  = (khash_t(s2i)*)lix->dict;

    int32_t v = 0;
    int32_t name_offset = 0;
    while (name_offset < length_of_names)
    {
        const char *name = names + name_offset;
        name_offset += strlen(name) + 1;
        k = kh_put(s2i, d, name, &absent);
        if (absent < 0) {
            hts_log_error("Out of memory!");
            goto fail;
        } else if (absent) {
            char *name_dup = strdup(name);
            if (name_dup) {
                kh_key(d, k) = name_dup;
                kh_val(d, k) = v++;
            } else {
                kh_del(s2i, d, k);
                goto fail;
            }
        }
    }

    lix->linear_offsets = (linear_offset_t *)calloc(
        lix->n, sizeof(linear_offset_t));
    if (!lix->linear_offsets) {
        hts_log_error("Out of memory!");
        goto fail;
    }

    int i;
    for (i = 0; i < lix->n; ++i) {
        linear_offset_t *loffset = &lix->linear_offsets[i];
        ret = check_bgzf_read(
            bgzf_read(lix_fp, &loffset->n, sizeof(loffset->n)),
            sizeof(loffset->n));
        if (ret < 0) goto fail;
        loffset->m = loffset->n;
        loffset->offset = (uint64_t *)calloc(loffset->n, sizeof(uint64_t));
        if (!loffset->offset) {
            hts_log_error("Out of memory!");
            goto fail;
        }
        ret = check_bgzf_read(
            bgzf_read(lix_fp, loffset->offset, loffset->n*sizeof(uint64_t)),
            loffset->n*sizeof(uint64_t));
        if (ret < 0) goto fail;
    }

    bgzf_close(lix_fp);
    free(names);
    
    return lix;
    
    fail:
        lix_destroy(lix);
        bgzf_close(lix_fp);
        if (names) free(names);
        return NULL;
}


int lix_equal(const lix_t *lix_a, const lix_t *lix_b) {
    if (lix_a == NULL && lix_b == NULL) return 0;
    if (lix_a == NULL && lix_b != NULL) return 1;
    if (lix_a != NULL && lix_b == NULL) return 1;
    if (lix_a->fmt != lix_b->fmt) return 1;
    if (lix_a->min_shift != lix_b->min_shift) return 1;
    if (lix_a->n != lix_b->n) return 1;
    int i, j;
    for (i = 0; i < lix_a->n; ++i) {
        linear_offset_t *loffset_a = &lix_a->linear_offsets[i];
        linear_offset_t *loffset_b = &lix_b->linear_offsets[i];
        if (loffset_a->n != loffset_b->n) return 1;
        for (j = 0; j < loffset_a->n; ++j) {
            if (loffset_a->offset[j] != loffset_b->offset[j]) return 1;
        }
    }

    char **names_a = (char **)malloc(lix_a->n*sizeof(char *));
    if (!names_a) {
        hts_log_error("Out of memeory!");
        return -1;
    }
    char **names_b = (char **)malloc(lix_b->n*sizeof(char *));
    if (!names_b) {
        hts_log_error("Out of memeory!");
        return -1;
    }


    int ret = lix_get_names(lix_a, (const char **)names_a);
    if (ret < 0) {
        hts_log_error("Failed to get names from lix_a!");
        free(names_a);
        free(names_b);
        return -1;
    }
    ret = lix_get_names(lix_b, (const char **)names_b);
    if (ret < 0) {
        hts_log_error("Failed to get names from lix_b!");
        free(names_a);
        free(names_b);
        return -1;
    }

    for (i = 0; i < lix_a->n; ++i) {
        if (strcmp(names_a[i], names_b[i])) return 0;
    }

    free(names_a);
    free(names_b);
    return 0;
} 


lix_iter_t *lix_iter_init(void) {
    lix_iter_t *iter = (lix_iter_t *)calloc(1, sizeof(lix_iter_t));
    if (!iter) {
        hts_log_error("Out of memory!");
        exit(1);
    }
    return iter;
}


void lix_iter_destroy(lix_iter_t *iter) {
    if (iter) {
        if (iter->contig) free(iter->contig);
    }
    free(iter);
}


lix_iter_t *lix_query(const lix_t *idx, const char *contig, int64_t start, int64_t end)
{
    lix_iter_t *iter = lix_iter_init();
    if (!iter) return NULL;
    iter->has_next = 0;
    iter->contig = strdup(contig);
    if (iter->contig == NULL) {
        hts_log_error("Out of memory!");
        exit(1);
    }
    if (start < 0) start = 0;
    if (end < 0) end = 0;
    if (start > end) {
        hts_log_error("Start position is greater than end position.");
        exit(1);
    }
    ldiv_t div = ldiv(start, 1<<idx->min_shift);
    if (div.rem != 0) {
        hts_log_warning("Start position is not multiple of %d(1<<%d), "
            "query may be very inefficient.", 1<<idx->min_shift, idx->min_shift);
    }
    iter->start = start;
    iter->end = end;
    int tid = lix_name2id(idx, iter->contig);
    if (tid < 0) {
        hts_log_error("Contig %s not fount in the index.", iter->contig);
        exit(1);
    }
    if (div.quot >= idx->linear_offsets[tid].n) {
        iter->offset = 0;
        iter->has_next = 0;
        return iter;
    }
    iter->offset = idx->linear_offsets[tid].offset[start/(1<<idx->min_shift)];
    iter->has_next = 1;
    return iter;
}


int lix_iter_next(lix_iter_t *iter, htsFile *fp, const bcf_hdr_t *hdr, bcf1_t *rec)
{
    if (fp->is_bgzf == 0) {
        hts_log_error("Only bgzf compressed files can be used with iterators.");
        return -2;
    }

    if (iter->has_next == 0) {
        return -1;
    }

    int ret = 0;
    if (iter->offset) {
        if (bgzf_seek(fp->fp.bgzf, iter->offset, SEEK_SET) != 0) {
            hts_log_error("Failed to seek to offset %lu.", iter->offset);
            return -2;
        }
        iter->offset = 0; // only seek once
    }

    const char *rec_contig = NULL;
    while (1) {
        ret = bcf_read1(fp, hdr, rec);
        // reached the end of the file, or on a read error
        if (ret != 0) {
            iter->has_next = 0;
            return ret;
        }
        
        rec_contig = bcf_hdr_id2name(hdr, rec->rid);
        if (strcmp(rec_contig, iter->contig) == 0) {
            // out of the query region
            if (rec->pos > iter->end) {
                iter->has_next = 0;
                return -1;
            }

            // skip the record before the query region
            if (rec->pos + rec->rlen <= iter->start) {
                continue;
            }

            // record is overlap with the query region
            return 0;
        } else {
            // reached the next contig
            iter->has_next = 0;
            return -1;
        }
    }

    return 0;
}

