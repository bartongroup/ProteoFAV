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
