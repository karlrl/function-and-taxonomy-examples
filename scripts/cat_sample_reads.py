from collections import defaultdict
from glob import glob
import os

import click


MAPPING_FILE_PATH = './metadata/SRS_SRR_ids_linked.txt'


def concatenate_reads(sample_id, srr_paths, direction, output_dir):
    out_filename = '{}_R{}.fastq.gz'.format(sample_id, direction)
    with open(os.path.join(output_dir, out_filename), 'wb') as out_file:
        for path in srr_paths:
            with open(path, 'rb') as srr_file:
                out_file.write(srr_file.read())


@click.command()
@click.argument('fastqs_dir', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path(exists=True))
def cat_sample_reads(fastqs_dir, output_dir):
    sample_map = defaultdict(list)

    with open(MAPPING_FILE_PATH) as f:
        lines = f.readlines()
        lines = lines[1:]
        for line in lines:
            srs, srr = line.split()
            sample_map[srs].append(srr)

    click.echo('Concatenating reads for {} samples:'.format(len(sample_map)))

    for sample_id, srr_ids in sample_map.items():
        forward_reads = []
        reverse_reads = []

        for srr_id in srr_ids:
            forward_reads += glob(fastqs_dir + '/{}_1.fastq.gz'.format(srr_id))
            reverse_reads += glob(fastqs_dir + '/{}_2.fastq.gz'.format(srr_id))

        assert len(forward_reads) == len(reverse_reads)
        if forward_reads:
            click.echo('{} forwards...'.format(sample_id))
            concatenate_reads(sample_id, forward_reads, '1', output_dir)
            click.echo('{} reverse...'.format(sample_id))
            concatenate_reads(sample_id, reverse_reads, '2', output_dir)
        else:
            click.secho(
                'No paired reads found for "{}"'.format(sample_id),
                fg='yellow',
            )


if __name__ == '__main__':
    cat_sample_reads()
