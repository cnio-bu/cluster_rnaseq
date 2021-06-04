## Common pipeline level functions

def is_single_end(sample):
    return pd.isnull(units.loc[sample, 'fq2']).all()


def is_multi_lane(sample):
    return not pd.isnull(units.loc[sample, 'lane']).all()
