
import config

import csv
from collections import defaultdict
import time

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


def write_multivariant_result_for_curly(multivar_result, passports):
    field_paths_for_curly = [('classification', 'rank1'),
                        ('classification', 'rank2'),
                        ('country',),
                        #'morpho_type',
                        #'sw_group', 'tmp', 'region', 'het',
                        #'percent_haplos_close_to_ref',
                        #'sw_category'
                        ]

    passports_by_cat = defaultdict(dict)
    fields_for_curly = set()
    for sample_id, passport in passports.items():
        for cat_path_in_passport in field_paths_for_curly:
            curly_cat = cat_path_in_passport[-1]

            passport_item = passport
            for key in cat_path_in_passport:
                passport_item = passport_item.get(key, {})

            if passport_item:
                if isinstance(passport_item, dict):
                    raise ValueError('passport item should not be a dict, but a str')

                value = passport_item
                fields_for_curly.add(curly_cat)
                passports_by_cat[curly_cat][sample_id] = value

    passports_for_curly = {cat: samples for cat, samples in passports_by_cat.items() if cat in fields_for_curly}

    multivar_result['projections'] = DataFrame(multivar_result['projections'],
                                               index=multivar_result['samples'])

    multivar_dir = config.MULTIVAR_DIR
    multivar_dir.mkdir(exist_ok=True)
    curly_path =  multivar_dir / 'pcoa.curly'
    write_pca_curlywhirly_file(multivar_result,
                               curly_path,
                               categories=passports_for_curly)

    back_dir = multivar_dir / 'back'
    back_dir.mkdir(exist_ok=True)
    datestamp = time.strftime("%Y-%m-%d_%H:%M:%S", time.gmtime())
    curly_path_back = back_dir / f'pcoa.{datestamp}.pcoa'
    write_pca_curlywhirly_file(multivar_result,
                               curly_path_back,
                               categories=passports_for_curly)


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
