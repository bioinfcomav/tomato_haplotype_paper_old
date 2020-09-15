import config

import gzip
import tempfile
import subprocess
from collections import defaultdict, OrderedDict

import numpy
import pandas

from variation import CHROM_FIELD, POS_FIELD


class Genes:
    def __init__(self):
        self._get_genes()

    def _get_genes(self):
        gene_ranges = []
        ids = []
        genes = {}
        for line in gzip.open(config.GENE_MODELS_GFF, 'rt'):
            if line.startswith('#'):
                continue

            items = line.split()
            if items[2] != 'gene':
                continue

            chrom = items[0]
            start = int(items[3])
            end = int(items[4])
            annotation = dict(item.split('=') for item in items[8].split(';'))
            id_ = annotation['ID'].split(':')[-1].split('.')[0]
            name = annotation['Name']
            genes[name] = {}
            gene_ranges.append((chrom, start, end))
            ids.append(id_)

        for line in gzip.open(config.GENE_MODELS_GFF, 'rt'):
            if line.startswith('#'):
                continue

            items = line.split('\t')
            if items[2] != 'mRNA':
                continue
            annotation = dict(item.split('=') for item in items[8].split(';'))
            parent_gene = annotation['Parent'].split(':')[-1].split('.')[0]

            try:
                parent_gene = genes[parent_gene]
            except KeyError:
                genes[parent_gene] = {}
                parent_gene = genes[parent_gene]
    
            mrna_id = annotation['ID'].split(':')[-1].split('.')[0]
            function = annotation.get('Note')

            if function == 'Unknown':
                function = None

            if 'mrnas' not in parent_gene:
                parent_gene['mrnas'] = {}
            parent_gene['mrnas'][mrna_id] = {}
            if function:
                parent_gene['mrnas'][mrna_id]['function'] = function

        genes_df = pandas.DataFrame(gene_ranges,
                                    columns=('Chromosome', 'Start', 'End'), index=ids)
        self._genes_df = genes_df

        self.gene_annotations = genes

    def get_gene(self, gene_id):
        gene = self._genes_df.loc[gene_id, :]
        gene = dict(gene.iteritems())
        return gene

    def get_annotated_function(self, gene_id):
        functions = [mrna_annotation['function'] for mrna_annotation in self.gene_annotations[gene_id].get('mrnas', {}).values() if mrna_annotation.get('function')]
        if functions:
            return functions[0]
        else:
            return None

    def _write_bed(self, df, sorted=False):
        df['End+1'] = df['End'] + 1
        df['name'] = df.index

        if not sorted:
            df = df.sort_values(by=['Chromosome', 'Start'])

        fhand = tempfile.NamedTemporaryFile(suffix='.bed', mode='wt')
        df.to_csv(fhand, sep='\t', header=False, columns=['Chromosome', 'Start', 'End+1', 'name'], index=False)
        fhand.flush()
        return fhand

    def get_closest_genes_to_vars(self, variations):
        df = {'Chromosome': numpy.array([chrom.decode() for chrom in variations[CHROM_FIELD]]),
              'Start': variations[POS_FIELD],
              'End': variations[POS_FIELD]}
        df = pandas.DataFrame(df)

        df_fhand = self._write_bed(df)
        genes_fhand = self._write_bed(self._genes_df)

        cmd = ['bedtools', 'closest', '-b', genes_fhand.name, '-a', df_fhand.name]
        process = subprocess.run(cmd, check=True, capture_output=True, text=True)
        genes = {line.split()[-1] for line in process.stdout.splitlines()}
        return genes

    def get_dist_to_closest_gene_for_vars(self, variations):
        df = {'Chromosome': numpy.array([chrom.decode() for chrom in variations[CHROM_FIELD]]),
              'Start': variations[POS_FIELD],
              'End': variations[POS_FIELD]}
        df = pandas.DataFrame(df)

        df_fhand = self._write_bed(df)
        genes_fhand = self._write_bed(self._genes_df)

        cmd = ['bedtools', 'closest', '-b', genes_fhand.name, '-a', df_fhand.name]
        process = subprocess.run(cmd, check=True, capture_output=True, text=True)

        dists = OrderedDict()
        for line in process.stdout.splitlines():
            items = line.split()
            var_start = int(items[1])
            var_end = int(items[2])
            gene_start = int(items[5])
            gene_end = int(items[6])
            chrom = items[0]

            dist_to_gene = min((abs(var_start - gene_start), abs(var_start - gene_end), abs(var_end - gene_start), abs(var_end - gene_end)))
            try:
                prev_dist = dists[(chrom, var_start)]
                if prev_dist > dist_to_gene:
                    dists[(chrom, var_start)] = dist_to_gene
            except KeyError:
                dists[(chrom, var_start)] = dist_to_gene

        return list(dists.values())

    @property
    def gene_ids(self):
        return iter(self._genes_df.index)

    def get_genes_in_region(self, chrom, start, end):
        df = {'Chromosome': [chrom],
              'Start': numpy.array([start], dtype=int),
              'End': numpy.array([end], dtype=int)}
        df = pandas.DataFrame(df)
        df.index = ['region']
        df_fhand = self._write_bed(df)

        genes_fhand = self._write_bed(self._genes_df)

        cmd = ['bedtools', 'intersect', '-b', genes_fhand.name, '-a', df_fhand.name, '-wb']
        process = subprocess.run(cmd, check=True, capture_output=True, text=True)

        genes = []
        for line in process.stdout.splitlines():
            gene_id = line.split('\t')[-1]
            genes.append(self.get_gene(gene_id))

        return genes


def read_go_annotation():
    first = True
    go_annotation = defaultdict(list)
    for line in config.GO_ANNOTATION.open():
        if first:
            first = False
            continue
        gene_id, go_term = line.split()
        gene_id = gene_id[:-4]
        go_annotation[gene_id].append(go_term)
    return go_annotation


def write_tomato_go_annotation_for_topgo():
    fhand = config.GO_ANNOTATION_TOPGO.open('wt')
    for gene, go_terms in read_go_annotation().items():
        go_terms = ', '.join(go_terms)
        fhand.write(f'{gene}\t{go_terms}\n')
    fhand.flush()


if __name__ == '__main__':
    genes = Genes()
    genes.get_genes_in_region('SL4.0ch02', 1, 1e5)
    genes.get_annotated_function('Solyc11g008780')

    #write_tomato_go_annotation_for_topgo()