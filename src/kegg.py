

import config

import requests
import pickle, gzip, hashlib

from pprint import pprint

from blast import get_cdna_ids_by_blasting


def get_pathway_info(path_id, organism, cache_dir=None):
    #print(f'Getting info for pathway {path_id}')

    content_lines = request_kegg_info(f'http://rest.kegg.jp/get/{path_id}',
                                      cache_dir=cache_dir)

    info_keys = {'MODULE': 'modules', 'GENE': 'genes', 'COMPOUND': 'compounds',
                 'REF_PATHWAY': 'related_pathways'}

    info = {val: {} for val in info_keys.values()}

    for line in content_lines:
        if line[0] != ' ':
            section = line.split(' ', 1)[0]
            line = line[len(section):]
        line = line.strip()
        if section in info_keys:
            items = line.split(' ', 1)
            if len(items) == 1:
                id_, description = items[0], None
            elif len(items) == 2:
                id_, description = items
            else:
                print(content_lines)
                print(section)
                print(line)
                raise RuntimeError()
            key = info_keys[section]
            if key == 'genes':
                id_ = f'{organism}:{id_}'

            info[key][id_] = description
    return info


def get_pathways(organism='sly', add_solgenomic_ids=False, cache_dir=None):

    content_lines = request_kegg_info(f'http://rest.kegg.jp/list/pathway/{organism}',
                                      cache_dir=cache_dir)

    pathways = {}
    genes = set()
    for line in content_lines:
        path_id, description = line.split('\t')
        pathways[path_id] = {'description': description}
        pathways[path_id]['info'] = get_pathway_info(path_id, organism,
                                                     cache_dir=cache_dir)
        genes_in_pathway = pathways[path_id]['info'].get('genes', {}).keys()
        genes.update(genes_in_pathway)

        if add_solgenomic_ids:
            pathways[path_id]['solgenomic_gene_ids'] = {gene['solgenomic_id'] for gene in get_genes(genes_in_pathway, cache_dir=cache_dir) if gene.get('solgenomic_id')}

    return {'pathways': pathways, 'genes': genes}


def request_kegg_info(url, cache_dir=None):

    if cache_dir:
        key = hashlib.md5(str(url).encode()).hexdigest()
        cache_path = cache_dir / ('kegg_rest_response' + key + '.pickle')
        if cache_path.exists():
            return pickle.load(gzip.open(cache_path, 'rb'))

    response = requests.get(url)
    content_lines = response.content.decode().splitlines()

    if cache_dir:
        pickle.dump(content_lines, gzip.open(cache_path, 'wb'))
    return content_lines


def get_gene(gene_id, cache_dir=None):
    content_lines = request_kegg_info(f'http://rest.kegg.jp/get/{gene_id}',
                                      cache_dir=cache_dir)

    res = {'kegg_gene_id': gene_id}
    section = None
    protein = ''
    dna = ''
    for line in content_lines:

        for key in ['NCBI-ProteinID', 'NCBI-GeneID', 'UniProt']:
            if f'{key}:' in line:
                res[key.lower()] = line.strip().split(':')[-1]
            if 'AASEQ' in line:
                section = 'protein_seq'
                continue
            if section == 'protein_seq':
                if line[0] != ' ':
                    section = None
                else:
                    protein += line.strip().strip()
            if 'NTSEQ' in line:
                section = 'nucleotide_seq'
                continue
            if section == 'nucleotide_seq':
                if line[0] != ' ':
                    section = None
                else:
                    dna += line.strip().strip()
    if protein:
        res['protein_seq'] = protein
    if dna:
        res['nucleotide_seq'] = dna
        solgenomic_cdna_ids = get_cdna_ids_by_blasting(dna)
        if solgenomic_cdna_ids:
            solgenomic_ids = {id_.split('.')[0] for id_ in solgenomic_cdna_ids}
            res['ids_for_similar_solgenomic_genes'] = solgenomic_ids
            if len(solgenomic_ids) == 1:
                res['solgenomic_id'] = list(solgenomic_ids)[0]
                #print(res['solgenomic_id'])
    return res


def get_genes(kegg_gene_ids, cache_dir=None):

    kegg_gene_ids = sorted(kegg_gene_ids)
    if cache_dir:
        key = str(kegg_gene_ids)
        key = hashlib.md5(str(key).encode()).hexdigest()
        cache_path = cache_dir / ('kegg_genes' + key + '.pickle')
        if cache_path.exists():
            return pickle.load(gzip.open(cache_path, 'rb'))

    genes = []
    for gene in kegg_gene_ids:
        gene = get_gene(gene, cache_dir=config.CACHE_DIR)
        genes.append(gene)

    if cache_dir:
        pickle.dump(genes, gzip.open(cache_path, 'wb'))
    return genes


if __name__ == '__main__':
    #pathways = get_pathway_info('path:sly00510')
    pathways = get_pathways(add_solgenomic_ids=True, cache_dir=config.CACHE_DIR)
    genes = pathways['genes']

    get_genes(genes, cache_dir=config.CACHE_DIR)

    #pprint(pathways['pathways'])
