
from collections import OrderedDict

LABELS = [
          ('slc_ma', 'SLC MA'),
          ('slc_mx_sinaloa', 'SLC Mx Sinaloa'),
          ('slc_mx', 'SLC Mx'),
          ('slc_ca', 'SLC CA'),

          ('sll_mx', 'SLL Mx'),
          ('sll', 'SLL'),

          ('slc_co', 'SLC Co'),

          ('slc_ec_n_600m', 'SLC Ec N 600m'),
          ('slc_ec_c_800m', 'SLC Ec C 800m'),
          ('slc_ec_s_1000m', 'SLC Ec N 1000m'),
          ('slc_ec_guayaquil', 'SLC Ec Guayaquil'),

          ('slc_pe_n_1000m', 'SLC Pe N 1000m'),
          ('slc_pe_n_400m', 'SLC Pe N 1000m'),
          ('slc_pe_c', 'SLC Pe C'),

          ('sp_ec_n_wet_forest', 'SP Ec N wet forest'),
          ('sp_ec_s_dry_forest', 'SP Ec S dry forest'),
          ('sp_pe_n_hills', 'SP Pe N hills'),
          ('sp_pe_n_inter-andean', 'SP Pe N Inter-Andean'),
          ('sp_pe_desert', 'SP Pe desert'),

          ('slc_world', 'SLC world'),
          ('sll_vint', 'SLL vintage'),

          ('sp_ec', 'SP Ec'),
          ('sp_pe', 'SP Pe'),
          ('sp_pe_inter-andean', 'SP Pe Inter-Andean'),
          ('sp_montane', 'SP Montane'),
          ('slc_ec', 'SLC Ec'),
          ('slc_pe', 'SLC Pe'),

          ('slc_pe_n', 'SLC Pe N'),
          ('slc_pe_s', 'SLC Pe S'),

          ('sp_pe_x_sp_ec', 'SP Pe x SP Ec'),
          ('sp_pe_desert_x_sp_pe_n_inter_andean', 'SP Pe x SP Pe'),
          ('sp_x_sl', 'SP x SL'),
          ('sp_x_sl_cherry_cultivars', 'SP x SL cherry cultivars'),

          ('sp_x_sp', 'SP x SP'),
          ('sll_old_cultivars', 'SLL old cultivars'),
          ('sll_modern', 'SLL modern'),
          ('sp_sl', 'SP-SL'),

          ('slc_small', 'SLC small'),
          ('slc_big', 'SLC big'),
          ('SLC small', 'SLC small'),
          ('SLC big', 'SLC big'),
          ('SP Montane', 'SP Montane'),

          (None, 'Unclassified'),
          ('None', 'Unclassified'),
          ('', 'Unclassified'),
          ]
LABELS = OrderedDict(LABELS)

def get_long_label(short_label):
    return LABELS.get(short_label, short_label)

HAPLO_LABELS = {'not_classified': 'Unclassified',
                'out_0': 'Outlier',
                'sl': 'hSL',
                'sp_ecu': 'hEc',
                'sp_peru': 'hPe',
                'outlier': 'Green introgression'
                }

DIVERSITY_INDEX_LABELS = {'haplo_diversity': 'Mean haplotype diversity'}