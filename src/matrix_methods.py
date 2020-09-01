
import config

import numpy

from variation.matrix.methods import iterate_matrix_chunks, vstack
from variation.matrix.stats import counts_and_allels_by_row
from variation import MISSING_INT


def unique(matrix, chunk_size=config.SNPS_PER_CHUNK):
    alleles = set()
    for chunk in iterate_matrix_chunks(matrix, chunk_size=chunk_size):
        alleles.update(numpy.unique(chunk))
    return sorted(alleles)


def counts_and_alleles_by_row(gts, missing_value=None, chunk_size=config.SNPS_PER_CHUNK, alleles=None):

    if alleles is None:
        alleles = unique(gts)

    if missing_value is not None:
        alleles = set(alleles).difference([missing_value])

    alleles = sorted(set(alleles))

    allele_counts = None
    chunks_allele_counts = (counts_and_allels_by_row(chunk, alleles=alleles)[0] for chunk in iterate_matrix_chunks(gts, chunk_size=chunk_size))
    allele_counts = vstack(chunks_allele_counts, missing_value=MISSING_INT)

    return allele_counts, alleles
