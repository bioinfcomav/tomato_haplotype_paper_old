
from collections import OrderedDict

LABELS = [
          ('slc_ma', 'SLC MA'),
          ('slc_mx_sinaloa', 'SLC MX Sinaloa'),
          ('slc_mx', 'SLC MX'),
          ('slc_ca', 'SLC CA'),

          ('sll_mx', 'SLL MX'),

          ('slc_co', 'SLC CO'),

          ('slc_ec_n_600m', 'SLC EC N 600m'),
          ('slc_ec_c_800m', 'SLC EC C 800m'),
          ('slc_ec_s_1000m', 'SLC EC N 1000m'),
          ('slc_ec_guayaquil', 'SLC EC Guayaquil'),

          ('slc_pe_n_1000m', 'SLC PE N 1000m'),
          ('slc_pe_n_400m', 'SLC PE N 1000m'),
          ('slc_pe_c', 'SLC PE C'),

          ('sp_ec_n_wet_forest', 'SP EC N wet forest'),
          ('sp_ec_s_dry_forest', 'SP EC S dry forest'),
          ('sp_pe_n_hills', 'SP PE N hills'),
          ('sp_pe_n_inter-andean', 'SP PE N Inter-Andean'),
          ('sp_pe_desert', 'SP PE desert'),

          ('slc_world', 'SLC world'),
          ('sll_vint', 'SLL vintage'),

          ('sp_ec', 'SP EC'),
          ('sp_pe', 'SP PE'),
          ('sp_pe_inter-andean', 'SP PE Inter-Andean'),
          ('slc_ec', 'SLC EC'),
          ('slc_pe', 'SLC PE'),

          ('sp_pe_x_sp_ec', 'SP PE x SP EC'),
          ('sp_pe_desert_x_sp_pe_n_inter_andean', 'SP PE x SP PE'),
          ('sp_x_sl', 'SP x SL'),
          ('sp_x_sl_cherry_cultivars', 'SP x SL cherry cultivars'),

          ('sp_x_sp', 'SP x SP'),
          ('sll_old_cultivars', 'SLL old cultivars'),
          ('sll_modern', 'SLL modern'),

          (None, 'Unclassified'),
          ('None', 'Unclassified'),
          ]
LABELS = OrderedDict(LABELS)

def get_long_label(short_label):
    return LABELS.get(short_label, short_label)

HAPLO_LABELS = {'not_classified': 'Unclassified',
                'out_0': 'Outlier',
                'sl': 'SL',
                'sp_ecu': 'EC',
                'sp_peru': 'SP PE',
                'outlier': 'Green introgression'
                }

DIVERSITY_INDEX_LABELS = {'haplo_diversity': 'Mean haplotype diversity'}