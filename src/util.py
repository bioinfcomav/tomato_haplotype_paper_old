
def dict_to_str(some_dict):
    keys = sorted(some_dict.keys(), key=str)
    dict_tuples = [(key, some_dict[key]) for key in keys]
    return str(dict_tuples)
