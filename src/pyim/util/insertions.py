

def subset_samples(insertions, samples, logger=None):
    warn = print if logger is None else logger.warning

    # Check for missing samples.
    ins_samples = set(insertions['sample'])
    for sample in samples:
        if sample not in ins_samples:
            warn('- Missing insertions for sample {}'.format(sample))

    # Actually subset insertions.
    return insertions.ix[ insertions['sample'].isin(set(samples))]
