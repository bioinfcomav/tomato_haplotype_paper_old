
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

    BASE_DIR = HOME_DIR / '/home/david/Disco/population_ultimate'
    MODULE_DIR = BASE_DIR / 'src'
    sys.path.insert(0, str(MODULE_DIR))

    DEVEL_DIR = HOME_DIR / '/home/david/Development'

    VARIATION_DIR = DEVEL_DIR / 'variation5'
    sys.path.insert(0, str(VARIATION_DIR))

SOURCE_DATA_DIR = BASE_DIR / 'source_data'

PASSPORTS_CSV = SOURCE_DATA_DIR / 'passports.csv'

print(BASE_DIR)
