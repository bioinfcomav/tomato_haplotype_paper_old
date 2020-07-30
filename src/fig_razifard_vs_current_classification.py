
import config

import pandas

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import passport
import pop_building
import plot
import labels


def _get_classifications(passports):
    samples = []
    rank1_classification = []
    razifard_classification = []
    for sample_name, passport in passports.items():
        rank1_pop = labels.LABELS[pop_building.get_classification_from_passport(config.RANK1, passport)]
        razifard_pop = pop_building.get_classification_from_passport(config.RAZIFARD, passport)
        if not razifard_pop:
            continue
        razifard_classification.append(razifard_pop)
        rank1_classification.append(rank1_pop)
        samples.append(sample_name)

    razifard_classification = pandas.Series(razifard_classification, index=samples)
    rank1_classification = pandas.Series(rank1_classification, index=samples)
    return {'razifard_classification': razifard_classification,
            'rank1_classification': rank1_classification}


if __name__ == '__main__':
    passports = passport.get_sample_passports()

    classifications = _get_classifications(passports)

    plot_path = config.FIG_RAZIFARD_VS_CURRENT_CLASSIFICATION

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)

    plot.plot_table_classification_comparison(classifications['razifard_classification'],
                                              classifications['rank1_classification'],
                                              axes=axes)

    fig.tight_layout()
    fig.savefig(plot_path)
