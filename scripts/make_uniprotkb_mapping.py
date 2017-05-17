from collections import defaultdict
import csv
import gzip


import click


MAPPING_TYPES = {
    'KEGG',
    'KO',
    'NCBI_TaxID',
    'UniProtKB-ID',
    'UniRef100',
    'UniRef50',
    'UniRef90',
}


@click.command()
@click.argument('uniprot_id_mapping', type=click.Path(exists=True))
@click.argument('ko_mapping', type=click.Path(exists=True))
@click.option('--output-mapping', type=click.Path(), default='mapping_file.tsv')
def make_uniprotkb_mapping(uniprot_id_mapping, ko_mapping, output_mapping):
    uniprot_id_mapping_file = gzip.open(uniprot_id_mapping, 'rt')
    id_mapping_reader = csv.reader(
        uniprot_id_mapping_file, delimiter='\t')
    output_map = defaultdict(lambda: defaultdict(list))
    ko_to_uniprotkb_ac = defaultdict(list)

    for i, row in enumerate(id_mapping_reader):
        uniprotkb_ac, mapping_type, mapped_id = row

        if mapping_type not in MAPPING_TYPES:
            continue
        elif mapping_type == 'KO':
            ko_to_uniprotkb_ac[mapped_id].append(uniprotkb_ac)

        output_row = output_map[uniprotkb_ac]
        output_row[mapping_type].append(mapped_id)

    uniprot_id_mapping_file.close()

    with open(ko_mapping) as ko_mapping_file:
        ko_mapping_reader = csv.reader(ko_mapping_file, delimiter='\t')
        next(ko_mapping_reader)  # Skip header

        for row in ko_mapping_reader:
            ko, pathways, modules = row
            uniprotkb_acs = ko_to_uniprotkb_ac[ko]
            for ac in uniprotkb_acs:
                if pathways != 'NA':
                    output_map[ac]['KEGG Pathways'] += pathways.split(',')
                if modules != 'NA':
                    output_map[ac]['KEGG Modules'] += modules.split(',')

    with open(output_mapping, 'wt') as f:
        output_writer = csv.writer(f, delimiter='\t')
        header = [
            'UniProtKB-AC',
            'KEGG',
            'KEGG Modules',
            'KEGG Pathways',
            'KO',
            'NCBI_TaxID',
            'UniProtKB-ID',
            'UniRef100',
            'UniRef50',
            'UniRef90',
        ]
        output_writer.writerow(header)
        for uniprotkb_ac, mappings in output_map.items():
            row = [uniprotkb_ac]
            row += [','.join(mappings.get(key)) for key in header[1:]],
            output_writer.writerow(row)

    click.echo('Generated mapping file.')


if __name__ == '__main__':
    make_uniprotkb_mapping()
