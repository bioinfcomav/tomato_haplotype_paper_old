
from scipy._lib.six import b
import config

import numpy

from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from collections import defaultdict

import morphological
import matplotlib_support
import labels

def _int_key(value):
    if value is None:
        return 1000000000
    else:
        return value


def plot_ordinal_trait_counts(axes, trait_counts, sorted_trait_forms, sorted_pops):


    x_values = numpy.arange(len(sorted_pops))
    width = x_values[1] - x_values[0]
    bottom = numpy.zeros((len(sorted_pops),))
    bars = []
    labels = []
    for trait_form in sorted_trait_forms:
        this_counts = numpy.array([trait_counts.get(pop, {}).get(trait_form, 0) for pop in sorted_pops])
        bar = axes.bar(x_values, this_counts, bottom=bottom, width=width, label=trait_form)
        bottom += this_counts
        bars.append(bar)
        labels.append(trait_form)
    return {'bars': bars, 'labels': labels}


def plot_morphological_characterization(ordinal_trait_counts, sorted_ordinal_traits=None,
                                        sorted_pops=None, represent_missing=False):

    trait_labels = morphological.get_trait_labels()

    if sorted_ordinal_traits is None:
        sorted_ordinal_traits = sorted(ordinal_trait_counts.keys())
    print(sorted_ordinal_traits)

    if sorted_pops is None:
        sorted_pops = sorted(set([pop for trait in sorted_ordinal_traits for pop in ordinal_trait_counts[trait].keys()]))

    fig = Figure((10, 20))
    FigureCanvas(fig) # Don't remove it or savefig will fail later

    gridspec = GridSpec(len(sorted_ordinal_traits), 2, figure=fig, width_ratios=[0.9, 0.1], wspace=0)

    for trait_idx, trait in enumerate(sorted_ordinal_traits):
        trait_counts = ordinal_trait_counts[trait]
        trait_forms = set([trait_form for pop in sorted_pops for trait_form in trait_counts[pop].keys()])
        if not represent_missing:
            trait_forms = trait_forms.difference([None])
        trait_forms = sorted(trait_forms, key=_int_key)

        axes = fig.add_subplot(gridspec[trait_idx, 0])
        res = plot_ordinal_trait_counts(axes, trait_counts, trait_forms, sorted_pops)

        legend_axes = fig.add_subplot(gridspec[trait_idx, 1])
        this_trait_labels = [trait_labels[trait][label] for label in res['labels']]
        legend_axes.legend(res['bars'], this_trait_labels, loc='upper left')
        matplotlib_support.set_axes_background(legend_axes)
        matplotlib_support.turn_off_both_axis(legend_axes)
        matplotlib_support.set_axis_color(legend_axes, color='white')

        axes.set_ylabel(f'{trait}\nFreq.')

        if trait_idx + 1 < len(sorted_ordinal_traits):
            matplotlib_support.turn_off_x_axis(axes)
        else:
            x_values = numpy.arange(len(sorted_pops))
            x_labels = [labels.LABELS[pop] for pop in sorted_pops]
            matplotlib_support.set_x_ticks(x_values, x_labels, axes, rotation=45)

    fig.savefig(config.FIG_MORPHOLOGICAL_TYPES_CHARACTERIZATION)


def calc_freqs_from_counts(counts, remove_missing=False):
    freqs = {}
    for trait, trait_counts in counts.items():
        freqs[trait] = {}
        for pop, trait_counts_for_pop in trait_counts.items():
            if remove_missing:
                trait_counts_for_pop = {trait_form: this_counts for trait_form, this_counts in trait_counts_for_pop.items() if trait_form is not None}

            tot_counts = sum(trait_counts_for_pop.values())
            freqs[trait][pop] = {trait_form: this_counts / tot_counts for trait_form, this_counts in trait_counts_for_pop.items()}
    return freqs


if __name__ == '__main__':
    
    morphotypes_for_accs = morphological.read_morphological_classification()
    accs_by_morphotype = defaultdict(list)
    for acc, morphotype in morphotypes_for_accs.items():
        accs_by_morphotype[morphotype].append(acc)

    counts = morphological.calc_ordinal_character_counts_per_acc_type(accs_by_morphotype)

    freqs = calc_freqs_from_counts(counts, remove_missing=True)

    sorted_pops = ['sp_pe', 'sp_intermediate', 'sp_ec', 'slc_ec', 'slc_small', 'slc_big', 'sll']
    sorted_traits = ['Fruit Size', 'Fruit Elon', 'F stripes', 'Rib', 'Peanut fruit',
                     'Style Exser', 'Style Curv', 'Petal pos', 'Petal width', 'Inflor len', 'Irreg inflor',
                     'Leaf margin', 'Leaf type',
                     'Stem hair', 'Stem width']
    plot_morphological_characterization(freqs, sorted_pops=sorted_pops, sorted_ordinal_traits=sorted_traits)
