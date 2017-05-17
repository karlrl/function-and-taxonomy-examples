from collections import defaultdict
import csv
import os
import re

import click


class Sample:
    def __init__(self, sample_id, counts_by_subject_id, num_alignments):
        self.sample_id = sample_id
        self.counts_by_subject_id = counts_by_subject_id
        self.num_alignments = num_alignments


class Alignment:
    def __init__(self, query_id, subject_id, percent_identity, e_value):
        self.query_id = query_id
        self.subject_id = subject_id
        self.percent_identity = percent_identity
        self.e_value = e_value

    @classmethod
    def parse_from_line(cls, line):
        alignment = line.split('\t')
        return Alignment(
            alignment[0],
            alignment[1],
            float(alignment[2]),
            float(alignment[10]),
        )


def aggregate_by_subject_id(rows):
    subject_counts = defaultdict(int)
    total_alignments = 0

    for row in rows:
        alignment = Alignment.parse_from_line(row)
        subject_counts[alignment.subject_id] += 1
        total_alignments += 1

    return total_alignments, subject_counts


def collect_subjects(samples):
    subjects = set()

    for sample in samples:
        for subject in sample.counts_by_subject_id.keys():
            subjects.add(subject)

    return subjects


def join_by_subject_id(samples):
    """
    Takes dict of dicts, ``relative_abundances`` and yields each relative
    abundance, joined by subject ID.
    """
    subjects_seen = sorted(list(collect_subjects(samples)))

    for subject in subjects_seen:
        row = [subject]
        for sample in samples:
            row.append(sample.counts_by_subject_id.get(subject, 0))

        yield row


def filter_samples_by_mapping_type(samples, mapping_file, from_type, to_type):
    for sample in samples:
        sample.filtered_counts_by_subject_id = defaultdict(int)
        sample.filtered_num_alignments = 0

    with open(mapping_file, 'rt') as f:
        header = next(f).strip().split('\t')
        from_index = header.index(from_type)
        to_index = header.index(to_type)

        for line in f:
            row = line.strip().split('\t')
            from_id = row[from_index]
            try:
                to_ids = [id_ for id_ in row[to_index].split(',') if id_ != '']
            except IndexError:
                continue

            if not (from_id and to_ids):
                continue
            for sample in samples:
                if from_id in sample.counts_by_subject_id:
                    count = sample.counts_by_subject_id[from_id]

                    # There are potentially many mappings for a given from_id.
                    # Naively assign the alignment counts to all of these
                    # mappings.
                    for to_id in to_ids:
                        sample.filtered_counts_by_subject_id[to_id] += count

                    # Use the true number of reads aligned for the total count
                    sample.filtered_num_alignments += count

    for sample in samples:
        sample.counts_by_subject_id = sample.filtered_counts_by_subject_id
        sample.num_alignments = sample.filtered_num_alignments
        del sample.filtered_counts_by_subject_id
        del sample.filtered_num_alignments

    return samples


def make_summary(samples, to_type=''):
    sorted_samples = sorted(samples, key=lambda s: s.sample_id)
    alignment_title = '{} Alignments'.format(to_type or '').lstrip()
    subjects_title = '{} Subjects'.format(to_type or '').lstrip()
    return [
        ['Summary'] + [s.sample_id for s in sorted_samples],
        [alignment_title] + [str(s.num_alignments) for s in samples],
        [subjects_title] + [str(len(s.counts_by_subject_id)) for s in samples]
    ]


@click.command()
@click.argument('input_files', nargs=-1, type=click.Path(exists=True))
@click.option('--output-file', type=click.Path(), default='abundances.spf')
@click.option('--summary-file', type=click.Path(), default='summary.spf')
@click.option('--mapping-file', type=click.Path())
@click.option('--from-type', type=click.STRING)
@click.option('--to-type', type=click.STRING)
@click.option('--cutoff', type=click.FLOAT, default=90.0)
def aggregate_alignments(
        input_files, output_file, mapping_file, from_type, to_type,
        summary_file, cutoff):

    if not input_files:
        click.secho('No input files could be found.', fg='yellow')
        exit(1)

    samples = []

    click.echo(
        'Extracting alignment counts from {} samples...'.format(
            len(input_files))
    )

    for path in input_files:
        sample_basename = os.path.basename(path)
        sample_id = re.search('[\w]*', sample_basename).group()

        with open(path, 'rt') as f:
            total_alignments, counts_by_subject_id = aggregate_by_subject_id(f)

        samples.append(
            Sample(sample_id, counts_by_subject_id, total_alignments)
        )

    if mapping_file:
        click.echo('Filtering alignment counts...')
        samples = filter_samples_by_mapping_type(
            samples, mapping_file, from_type, to_type)

    click.echo(
        'Writing alignment counts for all samples to {}...'
        .format(output_file)
    )

    with open(output_file, 'wt') as f:
        header = ['NAME'] + [s.sample_id for s in samples]
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        writer.writerows(join_by_subject_id(samples))

    with open(summary_file, 'wt') as f:
        summary_writer = csv.writer(f, delimiter='\t')
        summary_table = make_summary(samples, to_type)
        summary_writer.writerows(summary_table)


if __name__ == '__main__':
    aggregate_alignments()
