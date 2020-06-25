
import config

import tempfile
from collections import OrderedDict, Counter
import math
import itertools

import numpy

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import patches

from rpy2.robjects.packages import importr
ape = importr('ape')
pegas = importr('pegas')
rbase = importr('base')

import networkx as nx

from variation.variations import VariationsH5
from variation.gt_writers.fasta import write_fasta

import colors
import passport
import pop_building


class NotEnoughHaplotypesError(ValueError):
    pass


def calc_haplotype_network(variations, min_num_snp_for_network=20, remove_indels=True):

    if variations.num_variations < min_num_snp_for_network:
        raise ValueError('not enough SNPs')

    tmp_dir = config.TMP_DIR
    tmp_dir.mkdir(exist_ok=True)

    fasta_fhand = tempfile.NamedTemporaryFile(suffix='.fasta',
                                              delete=True,
                                              dir=tmp_dir)
    try_to_align_easy_indels = not(remove_indels)

    stats = write_fasta(variations, fasta_fhand, remove_indels=remove_indels,
                        try_to_align_easy_indels=try_to_align_easy_indels,
                        put_hyphens_in_indels=True,
                        remove_sites_all_N=True)
    fasta_fhand.flush()

    seqs_bin = ape.read_dna(fasta_fhand.name, format="fasta")
    haplotypes = pegas.haplotype(seqs_bin)

    haplo_composition = rbase.attr(haplotypes, 'index')
    n_items = rbase.length(haplo_composition)[0]
    if n_items < 2:
        raise NotEnoughHaplotypesError()

    network = pegas.haploNet(haplotypes)

    # conversion to python
    seq_names = list(rbase.labels(seqs_bin))
    haplo_composition = [list(haplo_composition[idx]) for idx in range(n_items)]

    haplo_labels1 = list(rbase.labels(haplotypes))
    haplo_labels2 = list(rbase.attr(network, 'labels'))
    assert haplo_labels1 == haplo_labels2

    sample_names = {idx + 1: name.split()[0][:-5] for idx, name in enumerate(seq_names)}

    network = numpy.array(network).astype(int)[:, :3]
    sample_sets = {}
    new_network = []
    for row_idx in range(network.shape[0]):
        node1, node2, distance = network[row_idx, :]
        seqs1 = haplo_composition[node1 - 1]
        seqs2 = haplo_composition[node2 - 1]
        samples1 = [sample_names[seq] for seq in seqs1]
        samples2 = [sample_names[seq] for seq in seqs2]

        sample_set1 = tuple(sorted(set(samples1)))
        if sample_set1 in sample_sets:
            sample_set1_idx = sample_sets[sample_set1]
        else:
            sample_set1_idx = len(sample_sets) + 1
            sample_sets[sample_set1] = sample_set1_idx

        sample_set2 = tuple(sorted(set(samples2)))
        if sample_set2 in sample_sets:
            sample_set2_idx = sample_sets[sample_set2]
        else:
            sample_set2_idx = len(sample_sets) + 1
            sample_sets[sample_set2] = sample_set2_idx

        new_network.append((sample_set1_idx, sample_set2_idx, {'distance': distance}))

    network = new_network
    sample_sets = {idx: samples for samples, idx in sample_sets.items()}

    graph = nx.Graph()
    graph.add_edges_from(network)

    return {'network_dict': network,
            'network_nx': graph,
            'sample_sets': sample_sets
           }


def find_clusters_in_network(network, k=3):
    clusters_iter = nx.algorithms.community.centrality.girvan_newman(network)
    clusters = list(itertools.islice(clusters_iter, k - 2, k - 1))[0]
    return clusters


def _count_pops_in_cluster(cluster, sample_sets, pop_for_samples):
    counts = Counter()
    for node_id in cluster:
        for sample in sample_sets[node_id]:
            pop = pop_for_samples[sample]
            counts[pop] += 1
    return counts


def _name_clusters(clusters, sample_sets, pop_for_samples):
    cluster_pop_counts = {}
    for idx, cluster in enumerate(clusters):
        counts = _count_pops_in_cluster(cluster, sample_sets, pop_for_samples)
        cluster_pop_counts[idx] = counts

    counts_for_sll_mex = [cluster_pop_counts[idx]['sll_mex'] for idx in range(len(cluster_pop_counts))]
    counts_for_sp_peru = [cluster_pop_counts[idx]['sp_peru'] for idx in range(len(cluster_pop_counts))]

    peruvian_cluster = numpy.argmax(numpy.array(counts_for_sp_peru))
    mexican_cluster = numpy.argmax(numpy.array(counts_for_sll_mex))

    cluster_names = {idx: None for idx in range(len(clusters))}
    if (not(isinstance(peruvian_cluster, numpy.ndarray) or isinstance(mexican_cluster, numpy.ndarray) and
        peruvian_cluster != mexican_cluster)):
        cluster_names[peruvian_cluster] = 'peruvian'
        cluster_names[mexican_cluster] = 'mexican'
        other_names_count = 1
        for idx in range(len(cluster_names)):
            if cluster_names[idx] is not None:
                continue
            name = f'other_{other_names_count}'
            cluster_names[idx] = name
            other_names_count += 1
    else:
        other_names_count = 1
        for idx in range(len(cluster_names)):
            name = f'unknown_{other_names_count}'
            cluster_names[idx] = name
            other_names_count += 1
    return {cluster_names[idx]: cluster for idx, cluster in enumerate(clusters)}


def plot_network(nx_network, sample_sets, axes, pop_for_samples, clusters=None,
                 edge_color_schema=None, pop_color_schema=None):

    if edge_color_schema is None:
        edge_color_schema = colors.ColorSchema()

    if pop_color_schema is None:
        pop_color_schema = colors.ColorSchema()

    circle_scale = 0.0001

    if clusters:
        node_cluster_belonging = {node_id: id_ for id_, cluster in clusters.items() for node_id in cluster}
    else:
        node_cluster_belonging = None

    node_pop_counts = OrderedDict()
    sizes = OrderedDict()
    for node in nx_network.nodes():
        node_pop_counts[node] = Counter([pop_for_samples.get(sample) for sample in sample_sets[node]])
        sizes[node] = sum(node_pop_counts[node].values())

    pos = nx.spring_layout(nx_network, weight='distance')
    #pos = nx.drawing.nx_agraph.graphviz_layout(nx_network)

    pos = {node_id: pos_ for node_id, pos_ in nx.spring_layout(nx_network).items()}

    #nx.draw(nx_network, pos=pos, with_labels=False, ax=axes, node_size=0)

    gray = '#aaaaaa'
    
    for edge in nx_network.edges:
        cluster0 = edge[0]
        cluster1 = edge[1]
        x1, y1 = pos[cluster0]
        x2, y2 = pos[cluster1]
        if node_cluster_belonging:
            cluster0_belonging = node_cluster_belonging[cluster0]
            cluster1_belonging = node_cluster_belonging[cluster1]
            if cluster0_belonging == cluster1_belonging:
                color = edge_color_schema[cluster0_belonging]
        else:
            color = gray
        axes.plot((x1, x2), (y1, y2), c=color, zorder=0)

    for node_id, (x_pos, y_pos) in pos.items():
        size = sizes[node_id] * circle_scale
        pops = node_pop_counts[node_id].keys()
        pie_colors = [pop_color_schema[pop] for pop in pops]
        sample_counts = [node_pop_counts[node_id][pop] for pop in pops]
        patches_, texts = axes.pie(sample_counts,
                                  radius=math.sqrt(size),
                                  center=(x_pos, y_pos),
                                  colors=pie_colors)
        for patch in patches_:
            patch.zorder = 1

    pops_present = sorted(set(pop_for_samples.values()), key=str)
    circles = []
    labels = []
    for pop in pops_present:
        color = pop_color_schema.colors_used[pop]
        if pop is None:
            pop = 'unclassified'
        circles.append(patches.Circle((0, 0), color=color))
        labels.append(pop)

    axes.legend(circles, labels, prop={'size': 6})



if __name__ == '__main__':

    min_num_snp_for_network = 14
    remove_indels = False

    vars_path = config.WORKING_PHASED_AND_IMPUTED_H5
    variations = VariationsH5(str(vars_path), 'r')

    variations = variations.get_chunk(slice(1000, 1100))

    res = calc_haplotype_network(variations,
                                 min_num_snp_for_network=min_num_snp_for_network,
                                 remove_indels=remove_indels)
    #print(res)

    passports = passport.get_sample_passports()
    pops_descriptions = {config.RANK1: config.ALL_POPS}
    pops = pop_building.get_pops(pops_descriptions, passports)
    pop_for_samples = {sample: pop for pop, samples in pops.items() for sample in samples}

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    clusters = find_clusters_in_network(res['network_nx'])
    clusters = _name_clusters(clusters, res['sample_sets'], pop_for_samples)

    plot_network(res['network_nx'], res['sample_sets'],
                 axes, pop_for_samples, clusters=None,
                 edge_color_schema=None,
                 pop_color_schema=colors.ColorSchema(colors.CLASSIFICATION_RANK1_COLORS))

    plot_path = config.TMP_DIR / 'network.svg'
    fig.tight_layout()
    fig.savefig(plot_path)
