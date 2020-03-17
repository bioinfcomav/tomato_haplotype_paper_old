
import csv

import numpy
from pandas import DataFrame

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from variation.variations.multivariate import do_pcoa


def write_curlywhirly_file(coordinates, path, categories=None):

    object_names = coordinates.index

    if categories is None:
        categories = {'cat': {name: 'object' for name in object_names}}

    cat_keys = list(categories.keys())

    projections = coordinates
    _, n_dims = projections.shape

    projection_axiss = ['axis{}'.format(n_dim)  for n_dim in range(1, n_dims + 1)]

    writer = csv.writer(path.open('wt'), delimiter='\t')

    cat_fields = ['categories:{}'.format(cat) for cat in cat_keys]
    fields = cat_fields + ['label'] + projection_axiss
    writer.writerow(fields)

    for object_idx, object_name in enumerate(object_names):
        object_point = projections.iloc[object_idx]
        object_cats = [categories[cat].get(object_name, 'unset') for cat in cat_keys]
        row = object_cats + [object_name] + list(object_point)
        writer.writerow(map(str, row))


def write_pca_curlywhirly_file(pca, curly_path, categories=None):
    write_curlywhirly_file(pca['projections'], curly_path, categories=categories)


def some_nan_in_numpy(array):
    return numpy.isnan(numpy.sum(array))


def do_pcoa_from_dists(dist_result, max_dims=5):

    multivar_result = do_pcoa(dist_result['dists'])

    multivar_result['projections'] = multivar_result['projections'][:, :max_dims]
    multivar_result['samples'] = dist_result['samples']

    return multivar_result


def do_pca(dframe, n_components=3):

    if some_nan_in_numpy(dframe.values):
        raise ValueError('nans are not allowed')

    std_dframe = StandardScaler().fit_transform(dframe)

    pca = PCA(n_components=n_components)

    projections = pca.fit_transform(std_dframe)
    
    n_digits = n_components // 10
    fstring = '{:0' + str(n_digits) + 'd}'
    pc_names = [fstring.format(idx) for idx in range(n_components)]
    

    projections = DataFrame(data=projections,
                            columns=pc_names,
                            index=dframe.index)
    explained_variance = [value * 100 for value in pca.explained_variance_ratio_]

    return {'projections': projections,
            'explained_variance': explained_variance}