
import config

import subprocess
from pathlib import Path


def prepare_blast_db(fasta_file, out_dir, db_type, skip_if_exists=False):
    fname = fasta_file.stem
    out_base_path = out_dir / fname

    if db_type == 'nucl':
        out_path = Path(str(out_base_path) + '.nhr')

    if skip_if_exists and out_path.exists():
        return

    cmd = ['makeblastdb', '-in', str(fasta_file), '-out', str(out_base_path),
           '-dbtype', db_type]
    subprocess.run(cmd, check=True)
    return {'db_path': out_base_path}


if __name__ == '__main__':

    blast_db_dir = config.CACHE_DIR / 'tomato_blast_db'

    prepare_blast_db(config.TOMATO_GENOME_FASTA, out_dir=blast_db_dir, db_type='nucl',
                     skip_if_exists=True)
