def add_prefix(values, prefix):
    return [prefix + v if not v.startswith(prefix) else v for v in values]


def remove_prefix(values, prefix):
    return [v[len(prefix):] if v.startswith(prefix) else v for v in values]
