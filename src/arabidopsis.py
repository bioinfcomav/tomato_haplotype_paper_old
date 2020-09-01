
import config

from Bio import SeqIO



def _read_loci():

    loci = {}
    gene_model_id_to_locus_id = {}
    for gene_model_id, seq in SeqIO.index(str(config.ARABIDOPSIS_CDNAS_FASTA), "fasta").items():
        locus_id = gene_model_id.split('.')[0]
        loci[locus_id] = {gene_model_id: {'cdna': str(seq.seq)}}
        gene_model_id_to_locus_id[gene_model_id] = locus_id

    for gene_model_id, pep in SeqIO.index(str(config.ARABIDOPSIS_PEPS_FASTA), "fasta").items():
        locus_id = gene_model_id_to_locus_id[gene_model_id]
        locus = loci[locus_id]
        try:
            gene_model = locus[gene_model_id]
        except KeyError:
            gene_model = {}
            locus[gene_model_id] = gene_model
        gene_model['pep'] = str(pep.seq)
    return loci, gene_model_id_to_locus_id


LOCI, GENE_MODEL_ID_TO_LOCUS_ID = _read_loci()


def get_locus(locus_id):
    return LOCI[locus_id]


def get_longest_pep(locus_id):
    locus = get_locus(locus_id)
    peps = list(reversed(sorted([gene_model['pep'] for gene_model in locus.values()],
                                key=lambda x: len(x))))
    return peps[0]


if __name__ == '__main__':

    print(get_longest_pep('AT3G24050'))