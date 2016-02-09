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


def parse_mutation(variation_table, column_name='residues'):
    """
    Parse AA mutation strings of the form [native]/[mutant1]/[mutant2]/...
    into to columns and combine with the supplied DataFrame.

    :param variation_table: A pandas.DataFrame with AA substitutions
    :param column_name: Name of column in variation_table that contains residue mutation
    :return: A pandas.DataFrame with the additional columns, 'from_aa' and 'to_aa'
    """
    residues = variation_table[column_name]
    from_aa = []
    to_aa = []
    for i in residues:
        parsed = i.split('/')
        from_aa.append(parsed[0])
        to_aa.append(parsed[1:])
    variation_table['from_aa'] = pd.Series(from_aa)
    variation_table['to_aa'] = pd.Series(to_aa)
    return variation_table