## Common pipeline level functions

def is_multi_lane(sample):
    return not pd.isnull(units.loc[sample, 'lane']).all()

def is_single_end_file():

    if pd.isnull(units.loc[:, 'fq2']).all() == True and pd.isnull(units.loc[:, 'fq2']).any() == True:
        return True
    
    elif pd.isnull(units.loc[:, 'fq2']).all() == False and pd.isnull(units.loc[:, 'fq2']).any() == False:
        return False
    
    elif pd.isnull(units.loc[:, 'fq2']).all() == False and pd.isnull(units.loc[:, 'fq2']).any() == True:
        sys.exit("\nThe data in units.tsv is not clearly Single-End or Paired-End. Please check the file units.tsv.\n") 

single_end = is_single_end_file()