
from os import remove
import config

import math
import random
from collections import defaultdict, Counter
from pprint import pprint
import csv

import pandas
import numpy

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import lat_lon
from pca import write_pca_curlywhirly_file
import passport
import colors
import labels
import matplotlib_support


ECUADOR_PROVINCES = 'Azuay, Bolívar, Cañar, Carchi, Chimborazo, Cotopaxi, El Oro, Esmeraldas, Galápagos, Guayas, Imbabura, Loja, Los Ríos, Manabí, Morona Santiago, Napo, Orellana, Pastaza, Pichincha, Santa Elena, Santo Domingo de Tsáchilas, Sucumbíos, Tungurahua, Zamora Chinchipe'
PERU_DEPARTMENTS = 'Amazonas, Áncash, Apurímac, Arequipa, Ayacucho, Cajamarca, Callao, Cusco, Huancavelica, Huánuco, Ica, Junín, La Libertad, Lambayeque, Lima, Loreto, Madre de Dios, Moquegua, Pasco, Piura, Puno, San Martín, Tacna, Tumbes, Ucayali'
ECUADOR_PROVINCES = [province.strip() for province in ECUADOR_PROVINCES.split(',')]
PERU_DEPARTMENTS = [department.strip() for department in PERU_DEPARTMENTS.split(',')]

COUNTRY_FOR_REGIONS = {province: 'ECU' for province in ECUADOR_PROVINCES}
COUNTRY_FOR_REGIONS.update({department: 'PER' for department in PERU_DEPARTMENTS})

MISSPELLINGS = {'Carchí': 'Carchi',
                'Manabi': 'Manabí',
                'Cuzco': 'Cusco',
                'Ancash': 'Áncash',
                'Zamora-Chinchipe': 'Zamora Chinchipe',
                'Morona-Santiago': 'Morona Santiago',
                'Apurimac': 'Apurímac',
                'Junin': 'Junín',
                'Huanuco': 'Huánuco',
               }

VALID_COUNTRIES = {'ECU', 'PER', 'MEX', 'COL', 'NIC', 'HND', 'CRI', 'SLV'}

TAXONS_FOR_EXCEL_SHEET = {0: 'SP',
                          1: 'SLC',
                          2: 'SLC',
                          3: 'SLC',
                          4: 'SLL'
                         }

TRAIT_COLUMNS = {'tipo hoja 1 pimpi; 2 intermedio; 3 no pimpi': 'leaf_type',
                 'borde foliolos 1 entero; 2 lobulado; 3 aserrado; 4 muy aserrado': 'leaflet_margin',
                 'pilosidad 1 nula; 2 intermedia; 3 alta': 'stem_hairiness',
                 'inflorescencia 1 muy larga; 2 larga; 3 intermedia; 4 corta; 5 irregular': 'inflorescence',
                 'posición pétalos 1 muy atrás; 2 entre 1 y 3; 3 medio; 4 no': 'petal_position',
                 'exerción estigmática 1 mucha; 2 media; 3 poca; 4 nula': 'style_exsertion',
                 'curvatura estilo 1 si; 2 no': 'style_curvature',
                 'anchura pétalos 1 anchos; 2 medios; 3 estrechos': 'petal_width',
                 'tamaño fruto 1 muy peque (Peru); 2 medio; 3 grande (norte ECU); 4 ceras 3 cm; 5 ceras grande': 'fruit_size',
                 'forma fruto 1 redondo; 2 cacah lig aplast; 3 lig aplast; 4 aplastado; 5 alargado': 'fruit_shape',
                 'rayas en el fruto 1 si; 2 no': 'stripy_fruit',
                 'anchura tallo 1 estrecho; 2 medio; 3 ancho': 'stem_width',
                 'acostillado 1 no; 2 leve; 3 moderado; 4 fuerte': 'ribbing',
                }

TRAIT_TYPES = {'leaf_type': {1: 'Pimpinellifolium', 2: 'intermediate', 3: 'standard'},
               'leaflet_margin': {1: 'entire', 2: 'lobulate', 3: 'serrate', 4: 'very serrated'},
               'stem_hairiness': {1: 'none', 2: 'intermidiate', 3: 'high'},
               'inflorescence': {1: 'very long', 2: 'long', 3: 'intermediate', 4: 'short', 5: 'irregular'},
               'inflorescence_ordinal' : {4: 'very long', 3: 'long', 2: 'intermediate', 1: 'short'},
               'petal_position': {1: 'folded back', 2: 'between back and medium', 3: 'medium', 4: 'none'},
               'style_exsertion': {1: 'high', 2: 'medium', 3: 'same level', 4: 'inserted'},
               'style_curvature': {1: 'yes', 2: 'no'},
               'petal_width': {1: 'wide', 2: 'medium', 3: 'thin'},
               'fruit_size': {1: 'sp small', 2: 'sp medium', 3: 'sp big', 4: 'slc', 5: 'bigger', 6: 'big', 7: 'biggest'},
               'fruit_shape': {1: 'round', 2: 'peanut', 3: 'slightly flattened', 4: 'flattened', 5: 'long'},
               'fruit_elongation': {1: 'round', 2: 'slightly flattened', 3: 'flattened', 4: 'long'},
               'stripy_fruit': {1: 'yes', 2: 'no'},
               'stem_width': {1: 'thin', 2: 'medium', 3: 'wide'},
               'ribbing': {1: 'none', 2: 'slight', 3: 'moderate', 4: 'strong'},
               'presence_of_irregular_inflorescence': {0: 'not present', 1: 'present'},
               'peanut_fruit_shape': {0: 'not present', 1: 'present'}
              }


ORDINAL_TRAITS = ['leaf_type', 'leaflet_margin', 'stem_hairiness', 'inflorescence_ordinal', 'petal_position',
                  'style_exsertion', 'style_curvature', 'petal_width', 'fruit_size', 'fruit_elongation',
                  'stripy_fruit', 'stem_width', 'ribbing', 'presence_of_irregular_inflorescence',
                  'peanut_fruit_shape']

TRAIT_ABREVIATIONS = {'leaf_type': 'Leaf type',
                      'leaflet_margin': 'Leaf margin',
                      'stem_hairiness': 'Stem hair',
                      'inflorescence_ordinal': 'Inflor len',
                      'petal_position': 'Petal pos',
                      'style_exsertion': 'Style Exser',
                      'style_curvature': 'Style Curv',
                      'petal_width': 'Petal width',
                      'fruit_size': 'Fruit Size',
                      'fruit_elongation': 'Fruit Elon',
                      'stripy_fruit': 'F stripes',
                      'stem_width': 'Stem width',
                      'ribbing': 'Rib',
                      'presence_of_irregular_inflorescence': 'Irreg inflor',
                      'peanut_fruit_shape': 'Peanut fruit'}

ORIGINAL_TRAITS = list(TRAIT_TYPES.values())

TRAITS_TO_REVERSE_SCALE = {'style_exsertion': {4: 1, 3: 2, 2: 3, 1: 4},
                           'style_curvature': {2: 1, 1: 2},
                           'stripy_fruit': {2: 1, 1: 2},
                           'petal_position': {1: 4, 2: 3, 3: 2, 4: 1},
                           'petal_width': {1: 3, 2: 2, 3: 1},
                           'inflorescence': {5: 0, 4: 1, 3: 2, 2: 3, 1: 4},
                          }


def is_nan(item):
    return isinstance(item, float) and math.isnan(item)


def _get_country(country, department):
    if is_nan(country):
        department = MISSPELLINGS.get(department, department)
        country = COUNTRY_FOR_REGIONS[department]
    assert country in VALID_COUNTRIES
    return country


def _get_acc_id(acc_id, bgv):
    if is_nan(bgv):
        assert acc_id
        return acc_id
    if isinstance(bgv, (float, int)):
        bgv = int(bgv)
        bgv = f'BGV{bgv:0>6}'
    acc_id = bgv
    return bgv


def _trait_value_to_number(value):
    if is_nan(value):
        return None

    if value == '1 y 2':
        int_value = random.choice([1, 2])
    elif value == '2 y 3':
        int_value = random.choice([2, 3])
    elif value == '4 y 5':
        int_value = random.choice([4, 5])
    elif value == '1 y 5':
        int_value = random.choice([1, 5])
    else:
        int_value = int(value)
    return int_value


def read_collecting_sources():

    reader = csv.DictReader(config.ECOSYSTEM_DATA.open('rt'), delimiter='\t')

    origins = {}
    for row in reader:
        accession_id = row['Accession'].lower()
        collecting_source = row['origen colecta']
        biological_status = row['tipo muestra']
        natural_checked = row['Natural comprobado: no afectado por humanos']
        ruderal_checked = row['ruderal comprobado: hábitat afectado por seres humanos']
        weedy_checked = row['Mala hierba comprobado']
        tolerated_checked = row['tolerada: surge expotáneamente (adventicia), pero la aprovechan (4, 1)']
        cultivated_in_garden = row['cultivada pequeño huerto']
        comercial_cultivation = row['cultivo comercial']
        market = row['market']

        if natural_checked:
            origin = 'natural'
        elif ruderal_checked or weedy_checked:
            origin = 'disturbed'
        elif weedy_checked and tolerated_checked:
            origin = 'weede_or_semi-cultivated'
        elif tolerated_checked:
            origin = 'semi-cultivated'
        elif cultivated_in_garden or comercial_cultivation:
            origin = 'cultivated'
        elif collecting_source == '21' and biological_status == '300':
            origin = 'cultivated'
        elif collecting_source == '22' and biological_status == '300':
            origin = 'cultivated'
        elif collecting_source == '23' and biological_status == '300':
            origin = 'cultivated'
        elif market:
            origin = 'market'
        elif collecting_source == '4' and biological_status == '1':
            origin = 'semi-cultivated'
        elif collecting_source == '1' and biological_status == '1':
            origin = 'natural_or_disturbed'
        elif collecting_source == '10' and biological_status == '100':
            origin = 'natural_or_disturbed'
        elif collecting_source == '1' and biological_status == '2':
            origin = 'disturbed'
        elif collecting_source == '4' and biological_status == '2':
            origin = 'semi-cultivated'
        elif collecting_source == '4' and biological_status == '4':
            origin = 'semi-cultivated_or_cultivated'
        elif collecting_source == '23' and biological_status == '200':
            origin = 'semi-cultivated_or_cultivated'
        elif collecting_source == '200' and biological_status == '23':
            origin = 'semi-cultivated_or_cultivated'
        elif collecting_source == '300' and biological_status == '21':
            origin = 'semi-cultivated_or_cultivated'
        elif collecting_source == '300' and biological_status == '23':
            origin = 'semi-cultivated_or_cultivated'
        elif collecting_source == '300' and biological_status == '22':
            origin = 'semi-cultivated_or_cultivated'
        elif collecting_source == '300':
            origin = 'semi-cultivated_or_cultivated'
        elif collecting_source == '2':
            origin = 'weed_or_semi-cultivated'
        elif biological_status == '4':
            origin = None
        elif collecting_source and biological_status:
            msg = f'Uknown collecting origin: {collecting_source}, {biological_status}'
            raise RuntimeError(msg)
        else:
            origin = None
        origins[accession_id] = origin

    return origins


def read_morphological_data():

    collecting_sources = read_collecting_sources()

    accs = {}
    for sheet_idx in range(5):
        sheet = pandas.read_excel(config.MORPHOLOGICAL_SOURCE_DATA, sheet_name=sheet_idx)

        for _, row in sheet.iterrows():
            acc = {}            
            acc['taxon'] = TAXONS_FOR_EXCEL_SHEET[sheet_idx]
            acc_id = _get_acc_id(row['Entrada'], row['BGV'])
            acc_id = acc_id.lower()
            acc['collection_id'] = row['Entrada']
            acc['collecting_source'] = collecting_sources[acc_id]

            location = {}
            location['country'] = _get_country(row['País'], row['Dpto'])
            lon = lat_lon.lon_to_deg(row['longitud'])
            if (not is_nan(lon)) and lon:
                location['longitude'] = lon
                location['latitude'] = lat_lon.lat_to_deg(row['latitud'])
            acc['location'] = location

            chracterization = {}            
            for trait_column, trait in TRAIT_COLUMNS.items():
                value = _trait_value_to_number(row[trait_column])
                value = TRAITS_TO_REVERSE_SCALE.get(trait, {}).get(value, value)
                chracterization[trait] = value
                if trait == 'inflorescence':
                    if value == 0:
                        chracterization['inflorescence_ordinal'] = None
                        chracterization['presence_of_irregular_inflorescence'] = 1
                    else:
                        chracterization['inflorescence_ordinal'] = value
                        chracterization['presence_of_irregular_inflorescence'] = 0
                if trait == 'fruit_shape':
                    if value in (1, 2):
                        if value == 2:
                            chracterization['peanut_fruit_shape'] = 1
                        else:
                            chracterization['peanut_fruit_shape'] = 0
                        value = 1
                    elif value is None:
                        value = None
                        chracterization['peanut_fruit_shape'] = None
                    else:
                        value = value - 1
                        chracterization['peanut_fruit_shape'] = 0
                    chracterization['fruit_elongation'] = value
            acc['characterization'] = chracterization
            accs[acc_id] = acc
    return accs


def get_morphological_table_for_ordinal_traits():
    data = defaultdict(list)
    acc_ids = []
    for acc_id, acc_data in read_morphological_data().items():
        acc_ids.append(acc_id)
        for trait in ORDINAL_TRAITS:
            data[trait].append(acc_data['characterization'][trait])
    morpho_data = pandas.DataFrame(data, index=acc_ids)
    return morpho_data


def fill_missing_data_with_means(dframe, max_missing_values=0):

    num_missing_values = dframe.isnull().sum(axis=1).values
    ok_rows = [row_id for row_id, num_missing in zip(dframe.index, num_missing_values) if num_missing <= max_missing_values]
    dframe = dframe.reindex(ok_rows)
    
    dframe.fillna(dframe.mean(), inplace=True)
    return dframe

    
def do_pca(dframe, num_dims=None):
    standarized_data = StandardScaler().fit_transform(dframe)
    pca = PCA(n_components=num_dims)
    principal_components = pca.fit_transform(standarized_data)
    result = {}
    old_dimensions = list(dframe.columns)
    new_dimensions = [f'dim_{idx + 1}' for idx in range(principal_components.shape[1])]
    result['projections'] = pandas.DataFrame(principal_components,
                                             index=list(dframe.index),
                                             columns=new_dimensions)
    result['var_percentages'] = pandas.Series(pca.explained_variance_ratio_ * 100,
                                              index=new_dimensions)
    result['princomps'] = pandas.DataFrame(pca.components_,
                                           index=new_dimensions,
                                           columns=old_dimensions)
    return result


def write_morphological_curlywirly(pca_result, morphological_data, passports, morpho_classification, morpho_classification_old=None):

    out_dir = config.MORPHOLOGICAL_DIR
    out_dir.mkdir(exist_ok=True)
    curly_path = out_dir / 'morpho_pca.curly'

    passports_for_curly = {}

    passports_for_curly['country'] = {acc_id: acc_data['location'].get('country') for acc_id, acc_data in morphological_data.items() if acc_data['location'].get('country')}
    passports_for_curly['taxon'] = {acc_id: acc_data['taxon'] for acc_id, acc_data in morphological_data.items()}

    passports_for_curly['rank1'] = {acc_id: passports.get(acc_id, {}).get('classification', {}).get('rank1') for acc_id in morphological_data.keys()}
    passports_for_curly['rank2'] = {acc_id: passports.get(acc_id, {}).get('classification', {}).get('rank2') for acc_id in morphological_data.keys()}
    passports_for_curly['geo_region'] = {acc_id: passports.get(acc_id, {}).get('geo_region') for acc_id in morphological_data.keys()}

    passports_for_curly['morpho'] = morpho_classification
    if morpho_classification_old:
        passports_for_curly['morpho_old'] = morpho_classification_old

    write_pca_curlywhirly_file(pca_result,
                               curly_path,
                               categories=passports_for_curly)


def read_morphological_classification(column='morpho_class'):
    return pandas.read_csv(config.CLASSIFICATIONS_MORPHOLOGICAL, sep=',', index_col=0, na_filter=False).to_dict()[column]


def write_morpho_csv(path, original_data, morpho_classification, write_collection_id=False):

    if write_collection_id:
        fields = ['Accession', 'CollectionID', 'Taxon', 'Country', 'Latitude', 'Longitude', 'Collecting_source']
    else:
        fields = ['Accession', 'Taxon', 'Country', 'Latitude', 'Longitude', 'Collecting_source']
    fields.append('Morphological classification')
    fields += sorted(TRAIT_TYPES.keys())
    writer = csv.DictWriter(path.open('wt'), fieldnames=fields)
    writer.writeheader()

    for acc_id, acc_data in original_data.items():
        row = {}
        row['Accession'] = acc_id.upper()
        row['Taxon'] = acc_data['taxon']
        location = acc_data.get('location')
        row['Country'] = location.get('country', '')
        row['Latitude'] = location.get('latitude', '')
        row['Longitude'] = location.get('longitude', '')
        if write_collection_id:
            row['CollectionID'] = acc_data['collection_id']
        characterization = acc_data.get('characterization', {})
        for trait in TRAIT_TYPES:
            value = characterization.get(trait, '')
            if value is None:
                value = ''
            row[trait] = str(value)
        row['Morphological classification'] = morpho_classification[acc_id]
        row['Collecting_source'] = acc_data['collecting_source']

        writer.writerow(row)


def plot_morphological_pca(pca_result, axes, classification, color_schema=None,
                           plot_principal_components=None,
                           principal_component_scaler=1):
    
    arrow_width = 0.01

    pops = defaultdict(list)
    for sample, klass in classification.items():
        pops[klass].append(sample)

    projections = pca_result['projections']
    artists = []
    labels = []
    for pop, samples in pops.items():
        samples_in_projections = set(samples).intersection(projections.index)
        this_pop_projections = projections.reindex(samples_in_projections)

        color = None if color_schema is None else color_schema[pop]
    
        artist = axes.scatter(this_pop_projections.values[:, 0],
                              this_pop_projections.values[:, 1], label=pop,
                              color=color)
        artists.append(artist)
        labels.append(pop)

    if 'var_percentages' in pca_result:
        axes.set_xlabel(f'Dim. 1 ({pca_result["var_percentages"][0]:.1f} %)')
        axes.set_ylabel(f'Dim. 2 ({pca_result["var_percentages"][1]:.1f} %)')

    if plot_principal_components is None:
        plot_principal_components = 'princomps' in pca_result
    
    if plot_principal_components:
        princomps = pca_result['princomps']

        vectors = princomps.T.iloc[:, :2]
        vector_xs = vectors.iloc[:, 0]
        vector_ys = vectors.iloc[:, 1]

        if principal_component_scaler != 1:
            modules = numpy.sqrt(vector_xs ** 2 + vector_ys ** 2)

            vector_xs = vector_xs * principal_component_scaler / modules
            vector_ys = vector_ys * principal_component_scaler / modules


        for trait, x_pos, y_pos in zip(princomps.columns, vector_xs, vector_ys):
            color = '#aaaaaa'
            axes.arrow(0, 0, x_pos, y_pos, color=color, width=arrow_width)
            axes.text(x_pos, y_pos, trait)
    return {'labels': labels, 'artists': artists}


def do_collecting_source_stats(morphological_data, remove_unknown_sources=True, calc_freqs=True, taxon_mapping=None):

    if taxon_mapping is None:
        taxon_mapping = {}
    counts = defaultdict(Counter)
    for acc_data in morphological_data.values():
        taxon = acc_data['taxon']
        collecting_source = acc_data['collecting_source']
        country = acc_data.get('location', {}).get('country')
        if remove_unknown_sources and (collecting_source is None or taxon is None or country is None):
            continue
        taxon_key = taxon_mapping.get((taxon, country), (taxon, country))
        #print(taxon_key)
        counts[taxon_key][collecting_source] += 1

    if calc_freqs:
        freqs = {}
        for taxon_country, this_counts in counts.items():
            tot_counts = sum(this_counts.values())
            freqs[taxon_country] = {source: source_counts / tot_counts for source, source_counts in this_counts.items()}
        counts = freqs

    return counts


def plot_collecting_sources(morphological_data, axes, color_schema, source_order, taxon_mapping, taxon_order):

    freqs = do_collecting_source_stats(morphological_data, taxon_mapping=taxon_mapping)

    sorted_taxons = sorted(freqs.keys(), key=lambda x: taxon_order.index(x))
    sources = {source for taxon_freqs in freqs.values() for source in taxon_freqs.keys()}
    sorted_sources = sorted(sources, key=lambda x: source_order.index(x))

    x_poss = numpy.arange(len(sorted_taxons))
    width = x_poss[1] - x_poss[0]

    bottoms = None
    artists = []
    labels = []
    for source in sorted_sources:
        heights = numpy.array([freqs[taxon].get(source, 0) for taxon in sorted_taxons])
        color = color_schema[source]
        artist =axes.bar(x_poss, heights, bottom=bottoms, width=width, label=source, color=color)
        artists.append(artist)
        labels.append(source)
        if bottoms is None:
            bottoms = heights
        else:
            bottoms += heights

    tick_labels = [f'{taxon} {country}' for taxon, country in sorted_taxons]
    matplotlib_support.set_x_ticks(x_poss, tick_labels, axes=axes, rotation=45)
    axes.set_ylabel('Freq.')
    return {'artists': artists, 'labels': labels}


def _to_int(value):
    if math.isnan(value) or value is None:
        return None
    else:
        return int(value)


def calc_ordinal_character_counts_per_acc_type(pops):

    morpho_data = get_morphological_table_for_ordinal_traits()
    morpho_data.columns = [TRAIT_ABREVIATIONS[trait] for trait in morpho_data.columns]

    accs_in_morpho_data = set(morpho_data.index)
    common_accs = {}
    for pop, accs in pops.items():
        common_accs[pop] = accs_in_morpho_data.intersection(accs)

    counts = {}
    for character, character_data in morpho_data.iteritems():
        counts[character] = {}
        for pop, accs in common_accs.items():
            this_pop_character_data = character_data.reindex(accs)
            counts[character][pop] = Counter([_to_int(value) for value in this_pop_character_data])

    return counts


def get_trait_labels():
    trait_labels = {}
    for trait, trait_forms in TRAIT_TYPES.items():
        trait_abreviation = TRAIT_ABREVIATIONS.get(trait, trait)

        if trait in TRAITS_TO_REVERSE_SCALE:
            trait_forms = {TRAITS_TO_REVERSE_SCALE.get(trait, {}).get(value, value): label for value, label in trait_forms.items()}
        
        trait_labels[trait] = trait_forms
        trait_labels[trait_abreviation] = trait_forms

    return trait_labels


if __name__ == '__main__':
    original_data = read_morphological_data()

    out_dir = config.MORPHOLOGICAL_DIR
    plot_path = out_dir / 'collecting_sources.svg'

    source_order = ['natural', 'natural_or_disturbed', 'disturbed',
                    'weed_or_semi-cultivated', 'semi-cultivated', 'semi-cultivated_or_cultivated',
                    'cultivated', 'market']

    taxon_mapping = {('SLC', 'CRI'): ('SLC', 'MA'),
                     ('SLC', 'MEX'): ('SLC', 'MA')}
    taxon_order = ('SP', 'PER'), ('SP', 'ECU'), ('SLC', 'COL'), ('SLC', 'CRI'), ('SLC', 'MEX'), ('SLC', 'MA'), ('SLC', 'PER'), ('SLC', 'ECU'), ('SLL', 'MEX')

    fig = Figure()
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)
    plot_collecting_sources(original_data, axes, color_schema=colors.ColorSchema(colors.SOURCE_COLORS),
                            source_order=source_order, taxon_mapping=taxon_mapping, taxon_order=taxon_order)
    matplotlib_support.set_axes_background(axes)
    axes.legend()
    fig.tight_layout()
    fig.savefig(plot_path)

    data = get_morphological_table_for_ordinal_traits()
    data.columns = [TRAIT_ABREVIATIONS[trait] for trait in data.columns]
    #print(Counter([acc_data['characterization']['stripy_fruit'] for acc_data in read_morphological_data().values()]))

    data = fill_missing_data_with_means(data, max_missing_values=8)
    pca_result = do_pca(data, num_dims=3)

    passports = passport.get_sample_passports()
    morpho_classification = read_morphological_classification()

    morpho_classification = {acc: klass for acc, klass in morpho_classification.items()}

    morpho_classification_old = read_morphological_classification('morpho_class_old')
    write_morphological_curlywirly(pca_result, original_data, passports, morpho_classification, morpho_classification_old)

    plot_path = out_dir / 'morphological_pca.svg'
    fig = Figure((10, 10))
    FigureCanvas(fig) # Don't remove it or savefig will fail later
    axes = fig.add_subplot(111)
    color_schema = colors.ColorSchema(colors.CLASSIFICATION_RANK1_COLORS)
    plot_morphological_pca(pca_result, axes, morpho_classification,
                           color_schema=color_schema,
                           plot_principal_components=True,
                           principal_component_scaler=3)
    matplotlib_support.set_axes_background(axes)
    axes.legend()
    fig.tight_layout()
    fig.savefig(str(plot_path))

    counts_per_pop = do_collecting_source_stats(original_data, remove_unknown_sources=True, calc_freqs=False, taxon_mapping=None)
    tot_num_accs = sum(count for counts in counts_per_pop.values() for count in counts.values())
    print('tot_num_accs', tot_num_accs)
    slc_ma_pops = [('SLC', 'CRI'), ('SLC', 'MEX')]
    slc_ec = [('SLC', 'ECU'), ('SLC', 'PER')]
    sp_pops = [('SP', 'ECU')]

    natural_or_disturbed_labels = ['natural_or_disturbed', 'natural', 'disturbed']
    cultivated = ['cultivated', 'semi-cultivated', 'market', 'semi-cultivated_or_cultivated', 'weed_or_semi-cultivated']

    pops = slc_ec
    labels = cultivated
    total_accs_in_pops_for_labels = 0
    total_accs_in_pops_for_labels_not = 0
    for pop in pops:
        total_accs_in_pops_for_labels += sum([count for label, count in counts_per_pop[pop].items() if label in labels])
        total_accs_in_pops_for_labels_not += sum([count for label, count in counts_per_pop[pop].items() if label not in labels])

    percent = total_accs_in_pops_for_labels / (total_accs_in_pops_for_labels + total_accs_in_pops_for_labels_not) * 100
    print(pops, labels, total_accs_in_pops_for_labels, total_accs_in_pops_for_labels_not, percent)