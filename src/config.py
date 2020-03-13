
import sys
from pathlib import Path
import getpass

HOME_DIR = Path.home()
USER = getpass.getuser()

if USER == 'jose':
    
    ANALYSES_DIR = HOME_DIR / 'magnet/analyses/'
    BASE_DIR = ANALYSES_DIR / 'ultimate_paper/'

    MODULE_DIR = BASE_DIR / 'src'
    sys.path.insert(0, str(MODULE_DIR))
    DEVEL_DIR = HOME_DIR / 'devel'
    VARIATION_DIR = DEVEL_DIR / 'variation5'
    sys.path.insert(0, str(VARIATION_DIR))
elif USER == 'david':

    BASE_DIR = HOME_DIR / '/home/david/Disco/ultimate_paper/'
    MODULE_DIR = BASE_DIR / 'src'
    sys.path.insert(0, str(MODULE_DIR))

    DEVEL_DIR = HOME_DIR / '/home/david/Development'

    VARIATION_DIR = DEVEL_DIR / 'variation5'
    sys.path.insert(0, str(VARIATION_DIR))

SOURCE_DATA_DIR = BASE_DIR / 'source_data'

PASSPORTS_CSV = SOURCE_DATA_DIR / 'passports.csv'

MAP_BACKGROUNDS_DIR = SOURCE_DATA_DIR / 'map_backgrounds'

NE_BACKGROUND_TIFF = MAP_BACKGROUNDS_DIR / 'NE1_50M_SR_W' / 'NE1_50M_SR_W.tif'
NE_BACKGROUND_CUT_PNG = MAP_BACKGROUNDS_DIR / 'ne1_50m_sr.cut.png'

FIGURES_DIR = BASE_DIR / 'figures'

GEOGRAPHIC_FIGURE_DIR = FIGURES_DIR / 'geographic_map'

HYPOTHESIS_PNG = GEOGRAPHIC_FIGURE_DIR / 'hypothesis.png'

CLASSIFICATIONS = BASE_DIR / 'classifications.csv'

CLASSIFICATION_RANKS = ['rank1', 'rank2', 'morpho_type', 'rank3']

SNPS_DIR = BASE_DIR / 'variations'

TIER1_PHASED_LOW_QUAL_09_MISSING_085 = SNPS_DIR / \
    '20190705.tier1.euchromatic.missing085_low_qual_09.phased.h5'
TIER1_PHASED_AND_IMPUTED_LOW_QUAL_09_MISSING_085 = SNPS_DIR / \
    '20190705.tier1.euchromatic.missing085_low_qual_09.phased_and_imputed.h5'
