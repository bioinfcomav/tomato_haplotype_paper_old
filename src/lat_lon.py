UNKOWN_CHAR = '-'
QUOTES = '"', "'"


def _latlon_to_deg(value):
    'It returns the latitude or longitude in degrees (float)'
    if isinstance(value, int):
        return float(value)
    elif isinstance(value, float) or value is None:
        return value
    elif isinstance(value, tuple) or isinstance(value, list):
        deg = float(value[0])
        minutes = float(value[1])
        sec = float(value[2])
        return deg + minutes / 60.0 + sec / 3600.0
    # So, it is a string
    value = value.strip()
    for quote in QUOTES:
        value = value.strip(quote)
    if not value or value[:-1] == (len(value) - 1) * UNKOWN_CHAR:
        return None

    def fix_unknown(string):
        'It replaces -- with 30'
        if string.startswith(UNKOWN_CHAR):
            return '30'
        elif string.endswith(UNKOWN_CHAR):
            return string[0] + '5'
        else:
            return string

    value = value.lower()
    if value[-1] in ('n', 's', 'w', 'e'):
        number = value[:-1]
        letter = value[-1]
        rev_digits = list(number)[::-1]
        rev_sec = rev_digits[:2]
        rev_min = rev_digits[2:4]
        rev_deg = rev_digits[4:]
        sec = int(fix_unknown(''.join(rev_sec[::-1])))
        minutes = int(fix_unknown(''.join(rev_min[::-1])))
        deg = int(''.join(rev_deg[::-1]))
        num = deg + minutes / 60.0 + sec / 3600.0
        if letter in ('s', 'w'):
            num = -num
        value = num
    else:
        value = float(value)
    return value


def lat_to_deg(value):
    'It returns a latitude in degrees'
#     if not value or set(value).difference(UNKOWN_CHAR):
#         raise ValueError('Empty latitude')

    deg = _latlon_to_deg(value)
    if deg and (deg < -90 or deg > 90):
        msg = 'Invalid latitude: %s -> %.2f' % (value, deg)
        raise ValueError(msg)
    return deg


def lon_to_deg(value):
    'It returns a longitude in degrees'
    deg = _latlon_to_deg(value)
    if deg and (deg < -180 or deg > 180):
        msg = 'Invalid longitude: %s -> %.2f' % (value, deg)
        raise ValueError(msg)
    return deg
