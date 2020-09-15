
import config

from tempfile import NamedTemporaryFile, TemporaryDirectory
from pathlib import Path
import subprocess
import gzip

from variation.gt_writers.vcf import write_vcf
from variation.gt_parsers.vcf import VCFParser
from variation.variations import VariationsArrays


BEAGLE_JAR = config.BEAGLE_JAR
BEAGLE_MEM = 10
DEFAULT_NE = 100000


def impute_variations(variations, genetic_map_path, beagle_mem=BEAGLE_MEM):

    with NamedTemporaryFile(mode='wb', suffix='.vcf') as vcf_fhand:
        write_vcf(variations, out_fhand=vcf_fhand)
        vcf_fhand.flush()
        imputed_vars = phase_and_impute_vcf_with_beagle(Path(vcf_fhand.name),
                                                        genetic_map_path=genetic_map_path,
                                                        beagle_mem=beagle_mem)
        vcf_fhand.close()
    return imputed_vars


def _sort_vcf(in_unsorted_vcf_path, out_sorted_vcf_path):
    sort_process = subprocess.Popen(['vcf-sort', in_unsorted_vcf_path],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = out_sorted_vcf_path.open('wb')
    gzip_process = subprocess.Popen(['bgzip'],
                                    stdin=sort_process.stdout,
                                    stdout=stdout,
                                    stderr=subprocess.PIPE)
    sort_process.stdout.close()
    gzip_process.communicate()
    assert gzip_process.returncode == 0
    stdout.close()


def phase_and_impute_vcf_with_beagle(vcf_path, genetic_map_path, imputed_vars=None,
                                     stdout=None, stderr=None, ap=True, gp=True,
                                     beagle_mem=BEAGLE_MEM):
    if stdout is None:
        stdout = subprocess.PIPE
    if stderr is None:
        stderr = subprocess.PIPE

    with TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        beagle_out = tmp_dir / 'beagle'
        cmd = ['java']
        cmd.append(f'-Xmx{beagle_mem}g')
        cmd.extend(['-jar', BEAGLE_JAR])
        cmd.append(f'map={genetic_map_path}')
        cmd.append(f'ne={DEFAULT_NE}')
        cmd.append(f'gt={vcf_path}')
        cmd.append(f'out={beagle_out}')
        cmd.append(f'ap={str(ap).lower()}')
        cmd.append(f'gp={str(gp).lower()}')

        beagle_out_vcf_path = Path(str(beagle_out) + '.vcf.gz')
        sorted_vcf = tmp_dir / 'beagle.sorted.vcf.gz'
        process = subprocess.run(cmd, stdout=stdout, stderr=stderr)
        if process.returncode:
            raise RuntimeError(process.stderr, process.stdout)

        _sort_vcf(beagle_out_vcf_path, sorted_vcf)
        vcf_parser = VCFParser(fhand=gzip.open(sorted_vcf))

        if imputed_vars is None:
            imputed_vars = VariationsArrays()

        imputed_vars.put_vars(vcf_parser)

        return imputed_vars
