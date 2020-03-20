
import numpy
import scipy
import pandas


def _align_pcas_using_procrustes(pcas, reference_projections=None):

    if reference_projections is None:
        reference_projections = pcas[0]['projections']

    for pca in pcas:
        projections = pca['projections']

        common_samples = reference_projections.index.intersection(projections.index)

        reference_projections_for_common_samples = reference_projections.loc[common_samples, :]
        projections_for_common_samples = projections.loc[common_samples, :]

        mtx1 = numpy.array(reference_projections_for_common_samples, dtype=numpy.double, copy=True)
        mtx1 -= numpy.mean(mtx1, 0)
        norm1 = numpy.linalg.norm(mtx1)
        if norm1 == 0:
           continue 
        mtx1 /= norm1

        mtx2 = numpy.array(projections_for_common_samples, dtype=numpy.double, copy=True)
        mtx2 -= numpy.mean(mtx2, 0)
        norm2 = numpy.linalg.norm(mtx2)
        mtx2 /= norm2

        mtx2_all = numpy.array(projections, dtype=numpy.double, copy=True)
        mtx2_all -= numpy.mean(mtx2_all, 0)
        norm2_all = numpy.linalg.norm(mtx2_all)
        if norm1 == 0:
           continue 
        mtx2_all /= norm2_all

        try:
            r_matrix, scale = scipy.linalg.orthogonal_procrustes(mtx1, mtx2)
        except ValueError:
            continue

        aligned_projections = numpy.dot(mtx2_all, r_matrix.T) * scale

        aligned_projections = pandas.DataFrame(aligned_projections,
                                               projections.index)
        yield {'projections': aligned_projections,
               'chrom': pca['chrom'],
               'win_start': pca['win_start']}


def _calc_pocas_mean_projections(pcoas):

    if not pcoas:
        raise ValueError('No PCoAs')

    if len(pcoas) == 1:
        return pcoas[0]['projections']

    projections_sum = None
    num_pcoas = 0
    for pcoa in pcoas:
        this_projections = pcoa['projections']

        if projections_sum is None:
            projections_sum = this_projections
        else:
            index = projections_sum.index.union(this_projections.index)

            projections_sum = projections_sum.reindex(index)
            this_projections = this_projections.reindex(index)

            projections_sum = projections_sum.add(this_projections, fill_value=0)
        num_pcoas += 1

    return pandas.DataFrame(projections_sum / num_pcoas, index=index)


def align_pcas_using_procrustes(pcoas):
    pcoas = list(pcoas)
    aligned_pcoas = _align_pcas_using_procrustes(pcoas)
    aligned_pcoas = list(aligned_pcoas)
    reference_projections = _calc_pocas_mean_projections(aligned_pcoas)
    del aligned_pcoas
    aligned_pcoas = _align_pcas_using_procrustes(pcoas, reference_projections)
    return aligned_pcoas
