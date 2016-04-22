import pandas as pd


def has_phenotype(variation_endpoint_json):
    results = pd.DataFrame()
    for variant, annotations in variation_endpoint_json.items():
        evidence = annotations['evidence']
        row = {'variant_id': variant}
        if 'Phenotype_or_Disease' in evidence:
            row.update({'has_phenotype_evidence': [True]})
        else:
            row.update({'has_phenotype_evidence': [False]})
        df = pd.DataFrame(row)
        results = results.append(df)
    return results


def variation_response_to_table(variation_response_json):
    """
    Convert the JSON response from the variation endpoint to a pandas.DataFrame

    :param variation_response_json: Decoded JSON from EnsEMBL variation query
    :return: pandas.DataFrame
    """
    results = pd.DataFrame()
    for variant, annotations in variation_response_json.items():

        parsed = dict([(k, [v]) for k, v in annotations.iteritems()])
        df = pd.DataFrame(parsed, index=[variant])
        results = results.append(df)
    results.index.name = 'variant_id'
    return results


def parse_mutation(variation_table, column_name='residues'):
    """
    Parse AA mutation strings of the form [native]/[mutant1]/[mutant2]/...
    into to columns and combine with the supplied DataFrame.

    :param variation_table: A pandas.DataFrame with AA substitutions
    :param column_name: Name of column in variation_table that contains residue mutation
    :return: A pandas.DataFrame with the additional columns, 'from_aa' and 'to_aa'
    """
    residues = variation_table[column_name].fillna('')
    from_aa = []
    to_aa = []
    for i in residues:
        parsed = i.split('/')
        from_aa.append(parsed[0])
        to_aa.append(parsed[1:])
    variation_table['from_aa'] = pd.Series(from_aa, index=residues.index)
    variation_table['to_aa'] = pd.Series(to_aa, index=residues.index)
    
    return variation_table


def parse_phenotypes_to_table(variation_table, drop_unparsed=False):
    """
    Parses the 'phenotypes' field from the EnsEMBL variation endpoint (a list of dictionaries) to a pandas.DataFrame.

    :param variation_table: A table with a 'phenotypes' column as produced by `variation_response_to_table`.
    :param drop_unparsed: If True, drop the original 'phenotypes' column from the result.
    :return: A pandas.DataFrame
    """

    # Create the phenotype table
    phenotype_table = pd.DataFrame()
    for variant, phenotypes in variation_table.phenotypes.iteritems():
        for d in phenotypes:
            parsed = dict([(k, [v]) for k, v in d.iteritems()])
            df = pd.DataFrame(parsed, index=[variant])
            phenotype_table = phenotype_table.append(df)
    phenotype_table.index.name = 'variant_id'

    # Merge the phenotype table back to the original input
    new_table = variation_table.merge(phenotype_table, how='left', left_index=True, right_index=True,
                                      suffixes=('', '_phenotypes'))

    if drop_unparsed:
        new_table.drop('phenotypes', axis=1, inplace=True)

    return new_table
