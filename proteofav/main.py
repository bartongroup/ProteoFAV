# -*- coding: utf-8

import sys
import click
import logging
import click_log

from proteofav.mergers import Tables

__all__ = ['main']

log = logging.getLogger('proteofav.config')


@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--pdb', type=str, default=None, help='Protein data bank identifier.')
@click.option('--chain', type=str, default=None, help='Protein structure chain identifier.')
@click.option('--uniprot', type=str, default=None, help='UniProt KnowledgeBase accession.')
@click.option('--output_type', default='csv',
              type=click.Choice(['csv', 'json', 'tab']),  # TODO add JALVIEW and chimera
              help='File format for the output.')
@click.option('-l', '--log', default=sys.stderr, help="Path to the logfile.",
              type=click.File('wb'))
@click.option('--add_dssp', is_flag=True,
              help="Whether to merge DSSP information to the output.")
@click.option('--add_annotations', is_flag=True,
              help="Whether to merge annotations information to the output.")
@click.option('--add_validation', is_flag=True,
              help="Whether to merge protein structure validation information to the output.")
@click.option('--add_variants', is_flag=True,
              help="Whether to merge genetic variant information to the output.")
# @click.option('--remove_redundant', is_flag=True,
#               help="Whether to remove columns with redundant information from the output.")
@click.argument('output', type=click.File('w'))
@click_log.simple_verbosity_option(default="INFO")
def main(pdb, chain, uniprot, output_type, log, add_dssp,
         add_annotations, add_validation, add_variants, output):
    """
    ProteFAV: a Python framework to process and integrate protein structure and
    features to genetic variants.
    OUTPUT: Path to the output file. Use `-` for stdout
    """

    logging.basicConfig(stream=log,
                        format='%(asctime)s - %(levelname)s - %(message)s ')

    table = Tables.generate(pdb_id=pdb,
                            uniprot_id=uniprot,
                            chains=chain,
                            sifts=True,
                            dssp=add_dssp,
                            validation=add_validation,
                            annotations=add_annotations,
                            variants=add_variants,
                            merge_tables=True,
                            sequence_check='ignore')

    if output_type == 'csv':
        table.to_csv(output, encoding='utf-8')
    elif output_type in {'jalview', 'chimera'}:
        raise NotImplementedError
    else:
        if hasattr(table, "to_" + output_type):
            fun = getattr(table, "to_" + output_type)
            fun(output)
    sys.exit(0)


@click.command()
@click.confirmation_option()
def setup():
    from proteofav.config import Defaults
    import os

    defaults = Defaults(os.path.join(os.path.dirname(__file__), 'config.ini'))
    default_db = os.path.expanduser("~/Downloads/")

    items = {k: v for k, v in defaults.__dict__.items() if k.startswith('db')}
    for attr, value in items.items():
        new_value = click.prompt(
            'Please enter a writable path for {}: '.format(attr),
            type=click.Path(writable=True),
            default=default_db)
        setattr(defaults, attr, new_value)

        if new_value != default_db:
            default_db = new_value

    new_email = click.prompt(
        'Please enter a valid email',
        type=str)
    if new_email:
        setattr(defaults, 'contact_email', new_email)

    defaults.write()
    log.info('Config file updated')
    sys.exit(0)


if __name__ == "__main__":
    pass
