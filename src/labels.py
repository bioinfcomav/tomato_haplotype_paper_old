
LABELS = {'sp_ec': 'SP EC',
          'sp_pe': 'SP PE',
          'sp_pe_inter-andean': 'SP PE Inter-Andean',
          'slc_ec': 'SLC EC',
          'sll_mx': 'SLL MX',
          'slc_ma': 'SLC MA',
          'slc_pe': 'SLC PE',
          'slc_co': 'SLC CO'
          }

def get_long_label(short_label):
    return LABELS.get(short_label, short_label)
