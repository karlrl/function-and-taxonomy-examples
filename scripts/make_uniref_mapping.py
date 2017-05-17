import csv
import gzip
from xml.etree import ElementTree as ET

import click


DB_FILENAME = 'uniref_mapping.db'
NAMESPACES = {'uniref': 'http://uniprot.org/uniref'}
ENTRY_TAG = '{{{}}}entry'.format(NAMESPACES['uniref'])


def parse_entry(entry):
    """Returns UniRef100, UniRef90 and UniRef50 IDs for a given ``entry``"""
    uniref100_id = entry.attrib['id']

    db_ref = entry.findall(
        './/uniref:representativeMember/uniref:dbReference',
        namespaces=NAMESPACES
    )
    db_ref = db_ref[0]

    uniref90_id = db_ref.findall(
        './/uniref:property[@type="UniRef90 ID"]', namespaces=NAMESPACES)
    if uniref90_id:
        uniref90_id = uniref90_id[0].attrib['value']
    else:
        uniref90_id = ''

    uniref50_id = db_ref.findall(
        './/uniref:property[@type="UniRef50 ID"]', namespaces=NAMESPACES)
    if uniref50_id:
        uniref50_id = uniref50_id[0].attrib['value']
    else:
        uniref50_id = ''

    return uniref100_id, uniref90_id, uniref50_id


@click.command()
@click.argument('uniref_xml_db', type=click.Path(exists=True, dir_okay=False))
@click.argument('output_file', type=click.Path(dir_okay=False))
def make_uniref_mapping(uniref_xml_db, output_file):
    row_count = 0
    progress_bar_batch_size = 10000
    progress_bar_kwargs = {
        'length': 106078676,  # Number of UniRef100 sequences (2017_04)
    }
    click.echo('Writing mapping file to {}'.format(output_file))

    with click.progressbar(**progress_bar_kwargs) as bar, \
            gzip.open(uniref_xml_db) as f_in, \
            open(output_file, 'wt') as f_out:

        parser = ET.iterparse(f_in, events=('start', 'end'))
        _, root = next(parser)
        output_writer = csv.writer(f_out, delimiter='\t')

        # Write the header for the output file
        output_writer.writerow(['UniRef100', 'UniRef90', 'UniRef50'])

        for event, element in parser:
            if element.tag == ENTRY_TAG and event == 'end':
                row = parse_entry(element)

                if len(row) < 3:
                    el = row[0] if len(row) > 0 else 'N/A'
                    click.secho(
                        'Element "{}" does not have mappings!'.format(el),
                        color='yellow',
                    )
                    continue

                row_count += 1
                output_writer.writerow(row)

                if row_count % progress_bar_batch_size == 0:
                    bar.update(progress_bar_batch_size)

                # Free up memory by removing the entry from the tree
                element.clear()
                root.remove(element)

        bar.update(row_count % progress_bar_batch_size)

    click.echo('Done. {} sequences mapped.'.format(row_count))


if __name__ == '__main__':
    make_uniref_mapping()
