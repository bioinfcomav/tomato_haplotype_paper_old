
import config


import numpy
import pandas


class NoStructureFiles(RuntimeError):
    pass


class IncompleteStructureRun(RuntimeError):
    pass


def _get_structure_path(dir_, suffix, base_fname=''):
    result_base_path = dir_ / base_fname
    structure_paths = [path for path in dir_.iterdir() if str(path).startswith(str(result_base_path))]
    paths = [path for path in structure_paths if path.suffix.endswith(suffix)]

    if len(paths) > 1:
        raise RuntimeError(f'Only one {suffix} file expected, but several found in: {dir_}')
    elif not paths:
        raise NoStructureFiles(f'The structure file was not found in: {dir_}')
    
    path = paths[0]
    return path


def _get_admixture_path(dir_):
    return _get_structure_path(dir_, '.meanQ', config.FASTSTRUCTURE_RESULT_BASE_FNAME)


def _get_log_path(dir_):
    return _get_structure_path(dir_, '.log', config.FASTSTRUCTURE_RESULT_BASE_FNAME)


def _parse_log_file(log_path):
    log_complete = False
    for line in log_path.open('rt'):
        if line.startswith('Total iterations = '):
            log_complete = True

        if line.startswith('Marginal Likelihood ='):
            marginal_likelihood = float(line.split('=')[1])

    if not log_complete:
        raise IncompleteStructureRun(f'There is an incomplete structure run in: {log_path}')

    return {'marginal_likelihood': marginal_likelihood}


def _parse_plink_fam_file(path):
    samples = [line.split()[0] for line in path.open('rt')]
    return {'samples': samples}


def _parse_admixture_file(path):
    admixtures = numpy.genfromtxt(path, dtype=float, delimiter='  ', skip_header=0)

    assert numpy.allclose(numpy.sum(admixtures, axis=1), 1)
    assert numpy.all(numpy.logical_not(numpy.isnan(admixtures)))

    return admixtures


def parse_faststructure_results():

    plink_fam_path = _get_structure_path(config.FASTSTRUCTURE_PLINK_DIR, '.fam')
    samples = _parse_plink_fam_file(plink_fam_path)['samples']

    ress = {}
    for path in config.FASTSTRUCTURE_RUN_DIR.iterdir():
        res = {}
        if not path.is_dir():
            continue

        try:
            log_path = _get_log_path(path)
        except (RuntimeError, NoStructureFiles):
            continue

        try:
            res['marginal_likelihood'] = _parse_log_file(log_path)['marginal_likelihood']
        except IncompleteStructureRun:
            continue

        try:
            admixture_path = _get_admixture_path(path)
        except (RuntimeError, NoStructureFiles):
            continue

        admixtures = pandas.DataFrame(_parse_admixture_file(admixture_path),
                                      index=samples)
        res['admixtures'] = admixtures

        k = admixtures.shape[1]
        res['k'] = k

        ress[k] = res
    return ress

if __name__ == '__main__':
    parse_faststructure_results()