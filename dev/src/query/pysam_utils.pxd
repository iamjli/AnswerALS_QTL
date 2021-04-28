from pysam.libchtslib cimport *


cdef extern from "htslib_util.h":

    # add *nbytes* into the variable length data of *src* at *pos*
    bam1_t * pysam_bam_update(bam1_t * b,
           	         size_t nbytes_old,
                         size_t nbytes_new,
                         uint8_t * pos)

    # now: static
    int aux_type2size(int)

    uint32_t * pysam_bam_get_cigar(bam1_t * b)
    uint32_t pysam_get_n_cigar(bam1_t * b)


from pysam.libcalignmentfile cimport AlignmentFile, AlignmentHeader, AlignmentHeader
from pysam.libcalignedsegment cimport makeAlignedSegment, AlignedSegment

# factory methods
cdef AlignedSegment makeAlignedSegment(
    bam1_t * src,
    AlignmentHeader header)

cdef inline uint32_t get_alignment_length(bam1_t *src):
    cdef uint32_t k = 0
    cdef uint32_t l = 0
    if src == NULL:
        return 0
    cdef uint32_t * cigar_p = bam_get_cigar(src)
    if cigar_p == NULL:
        return 0
    cdef int op
    cdef uint32_t n = pysam_get_n_cigar(src)
    for k from 0 <= k < n:
        op = cigar_p[k] & BAM_CIGAR_MASK
        if op == BAM_CSOFT_CLIP or op == BAM_CHARD_CLIP:
            continue
        l += cigar_p[k] >> BAM_CIGAR_SHIFT
    return l




