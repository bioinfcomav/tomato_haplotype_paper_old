
import config

import csv

import arabidopsis
import blast
import genome


# from Light Perception: A Matter of Time https://doi.org/10.1016/j.molp.2020.02.006

ARABIDOPSIS_CIRCADIAN_LIGHT_SIGNALING_MUTANTS = {
'hyh': 'AT3G17609',
'fhy3': 'AT3G22170',
'far1': ' AT4G15090',
'hy1': 'AT2G26670',
'cor27': 'AT5G42900',
'cor28': 'AT4G33980',
'hy5': 'AT5G11260',
'det1': 'AT4G10180',
'cop1': 'AT2G32950',
'pif5': 'AT3G59060',
'pif4': 'AT2G43010',
'pif3': 'AT1G09530',
'pif1': 'AT2G20180',
'uvr8': 'AT5G63860',
'ztl': 'AT5G57360',
'cry2': 'AT1G04400',
'cry1': 'AT4G08920',
'phye': 'AT4G18130',
'phyd': 'AT4G16250',
'phyc': 'AT5G35840',
'phyb': 'AT2G18790',
'phya': 'AT1G09570'}


def get_arabidopsis_circadian_annotated_loci():
    loci = []
    first = True
    for line in csv.reader(config.ARABIDOPSIS_CIRCADIAN_ANNOTATED_LOCI.open('rt'), delimiter='\t'):
        if first:
            first = False
            continue
        loci.append(line[0].split('/')[0])
    return loci


def _get_tomato_genes_similar_loci(arabiopsis_loci,
                                   cache_dir=config.CACHE_DIR,
                                   evalue_threshold=1e-40,
                                   only_first=False):

    genes = []
    genes_by_loci = {}
    arabidopsis_loci_by_tomato_gene = {}
    for locus_id in arabiopsis_loci:
        try:
            pep = arabidopsis.get_longest_pep(locus_id)
        except KeyError:
            continue
        similar_genes = blast.get_cdna_ids_by_blasting(pep, blast_program='tblastn',
                                                       evalue_threshold=evalue_threshold,
                                                       cache_dir=cache_dir)
        if only_first:
            if similar_genes:
                genes.append(similar_genes[0])
                genes_by_loci[locus_id] = similar_genes[0]
                arabidopsis_loci_by_tomato_gene[similar_genes[0]] = locus_id
        else:
            genes.extend(similar_genes)
            genes_by_loci[locus_id] = similar_genes
            for similar_gene in similar_genes:
                arabidopsis_loci_by_tomato_gene[similar_gene] = locus_id
    return {'genes': genes, 'genes_by_loci': genes_by_loci,
            'arabidopsis_loci_by_tomato_gene': arabidopsis_loci_by_tomato_gene}



def get_tomato_genes_similar_to_arabidopsis_circadian_annotated_loci(cache_dir=config.CACHE_DIR,
                                                                     evalue_threshold=1e-40,
                                                                     only_first=False):
    return _get_tomato_genes_similar_loci(get_arabidopsis_circadian_annotated_loci(),
                                          cache_dir=cache_dir,
                                          evalue_threshold=evalue_threshold,
                                          only_first=only_first)


def get_tomato_genes_similar_to_arabidopsis_circadian_mutants(cache_dir=config.CACHE_DIR,
                                                              evalue_threshold=1e-40,
                                                              only_first=False):

    return _get_tomato_genes_similar_loci(ARABIDOPSIS_CIRCADIAN_LIGHT_SIGNALING_MUTANTS.values(),
                                          cache_dir=cache_dir,
                                          evalue_threshold=evalue_threshold,
                                          only_first=only_first)


def _get_circadian_loci(cache_dir=config.CACHE_DIR, evalue_threshold=1e-30,
                       only_first=False, only_references=None):

    genes = genome.Genes()

    tomato_genes = {}
    for row in csv.DictReader(config.CIRCADIAN_LOCI.open('rt'), delimiter='\t'):
        arabidopsis_id = row['Arabidopsis id'].upper()

        if only_references:
            references = row['reference'].split(';')
            if not set(references).intersection(only_references):
                continue

        pep = None
        cdna = None
        if arabidopsis_id:
            pep = arabidopsis.get_longest_pep(arabidopsis_id)
        elif row['pep sequence']:
            pep = row['pep sequence']
        else:
            cdna = row['cDNA sequence']

        if pep:
            tomato_similar_genes = blast.get_cdna_ids_by_blasting(pep, blast_program='tblastn',
                                                                evalue_threshold=evalue_threshold,
                                                                cache_dir=cache_dir)
        else:
            tomato_similar_genes = blast.get_cdna_ids_by_blasting(cdna, blast_program='blastn',
                                                                  evalue_threshold=evalue_threshold,
                                                                  cache_dir=cache_dir)
        if only_first and tomato_similar_genes:
            tomato_similar_genes = [tomato_similar_genes[0]]

        for tomato_gene in tomato_similar_genes:
            gene = {'locus_id': tomato_gene.split('.')[0]}
            if arabidopsis_id:
                gene['similar_to_arabidopsis'] = arabidopsis_id
            #print(gene)
            tomato_genes[gene['locus_id']] = gene
    return tomato_genes


def get_tomato_genes_annotations():
    reader = csv.DictReader(config.TOMATO_GENE_ANNOTATIONS.open('rt'), delimiter='\t')
    annotations = {}
    for row in reader:
        annotation = {}
        annotation['molecular_functions'] = row['molecular_functions'].split(';')
        annotations[row['Gene id']] = annotation
    return annotations


def get_circadian_loci(cache_dir=config.CACHE_DIR, evalue_threshold=1e-30,
                       only_first=False, only_references=None):
    genes = _get_circadian_loci(cache_dir=cache_dir, evalue_threshold=evalue_threshold,
                               only_first=only_first, only_references=only_references)
    annotations = get_tomato_genes_annotations()
    for gene_id, gene_data in genes.items():
        try:
            annotation = annotations[gene_id]
        except:
            continue
        gene_data.update(annotation)
    return genes


if __name__ == '__main__':

    tomato_genes = get_circadian_loci(only_first=False, only_references=None)
    print(len(tomato_genes))

