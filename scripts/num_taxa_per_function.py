from __future__ import print_function, division
import argparse
import sys
from ete2 import NCBITaxa
ncbi = NCBITaxa()


def ncbi_taxa_lineage(tax_id, level2keep):

    '''Given NCBI taxon ID figure out ids for level of interest along
    lineage.'''

    # Get full taxonomic lineage based on NCBI taxon id.
    tax_lineage = ncbi.get_lineage(tax_id)
    tax_lineage_rank = ncbi.get_rank(tax_lineage)

    # Dictionary with levels as keys and NCBI ids as values.
    lineage = {}

    # Loop over keys/values of lineage rank dictionary and keep ids that
    # are at lineages of interest.
    for lineage_id, level in tax_lineage_rank.iteritems():
        if level == level2keep:

            if level in lineage:
                print("level " + level + " is already in lineage for "
                      + tax_id, file=sys.stderr)

            lineage[level] = lineage_id

    return(lineage)


def add_lineage_to_func(d, lin, uniref100_id, uniref50_id, uniref90_id, module,
                        pathway, ko, level2keep):

    '''Add taxon id under function and level of interest.'''

    cat = ["uniref100", "uniref90", "uniref50", "ko", "module",
           "pathway"]

    function_ids = [uniref100_id, uniref90_id, uniref50_id, ko, module,
                    pathway]

    # Intitialize index for "cat" list.
    cat_i = -1

    # Loop through functional ids.
    for full_func in function_ids:

        cat_i += 1

        # Skip if no function associated.
        if(full_func == ""):
            continue

        # Split on "," since sometimes multiple functions can be linked to
        # the same gene.
        for func in full_func.split(","):

            # Initialize empty dictionary for this function if undefined.
            if(func not in d[cat[cat_i]]):
                d[cat[cat_i]][func] = {}

            # Add lineage under function of interest.
            d[cat[cat_i]][func][lin[level2keep]] = None


def get_taxa_counts(d):

    '''Get counts of unique taxa for each level in order to output distribution
    of specificity of function within a given category.'''

    cat = ["uniref100", "uniref90", "uniref50", "ko", "module",
           "pathway"]

    # Dictionary containing all counts.
    counts = {}

    # Initialize the highest number of unique taxa observed thus far.
    # This is done so that the number of different bins to output is known
    # in advance.
    level_max = 0

    # Loop over all functional types.
    for c in cat:

        # Initialize empty dictionary for this type.
        counts[c] = {}

        # Loop over all individual functional ids under this functional type.
        for func in d[c].keys():

            # Get number unique taxa ids at level of interest.
            num_unique = len(d[c][func].keys())

            # Initialize this observed number of counts for this functional
            # type.
            if(num_unique not in counts[c]):
                counts[c][num_unique] = 0

            # Add an observation to this count of taxa for a functional id
            # under this functional type.
            counts[c][num_unique] += 1

            # Keep track of the ongoing highest number of unique taxa
            # for this level.
            if num_unique > level_max:
                level_max = num_unique

    return(counts, level_max)


def output_num_taxa(counts, level_max, level2keep):
    '''Output table of number of taxa that contain a given function for every
    taxonomic level.'''

    cat = ["uniref100", "uniref90", "uniref50", "ko", "module",
           "pathway"]

    outfile = open(level2keep+"_function_counts.txt", "w")

    # Print out headerline, with number of columns equal to the maximum
    # taxa counts observed.
    headerline = ["category"] + map(str, range(level_max+1))
    print("\t".join(headerline), file=outfile)

    # Loop over all functional types and print the distribution of taxa
    # counts for that type on 1 row.
    for c in cat:

        output_line = [c]

        for i in range(level_max+1):
            if i in counts[c]:
                output_line.append(str(counts[c][i]))
            else:
                # No function under this type is shared across this number of
                # taxa so set to 0.
                output_line.append(str(0))

        # Print out row.
        print("\t".join(output_line), file=outfile)

    outfile.close()


def main():

    parser = argparse.ArgumentParser(description="Parse mapping file to get\
             the taxonomic specificity of UniRef clusters and KEGG functions\
             at a specified taxonomic level. UniRef clusters at 50, 90,\
             and 100 percent are analyzed. KEGG orthologs, modules, and\
             pathways are also analyzed.")

    parser.add_argument("mapping_file", metavar="mapping_file", type=str,
                        nargs=1, help="File with mappings from UniRef clusters\
                        to NCBI taxa and KEGG levels.")

    parser.add_argument("level", metavar="level", type=str,
                        nargs=1, help="Level to get number of taxa that have\
                        function. One of superkingdom, kingdom, phylum,\
                        class, order, family, genus, species.")

    args = parser.parse_args()

    level2keep = args.level[0]

    # Read through mapping file and pull out 7 columns of interest:
    # UniRef 50, 90, and 100, KEGG orthologs, pathways, and modules, and NCBI
    # TaxIDs. Keep track of the number of different taxa that contains each
    # function.

    # First get dictionary which will contain keys for each type of function,
    # with each individual function id as a value. Each of these ids will
    # be the key to a deeper dictionary which will contain a set of different
    # taxa at each taxonomic level of interest.
    func2taxa = {}

    # Intitialize inner dictionaries for each function type.
    func2taxa["uniref100"] = {}
    func2taxa["uniref90"] = {}
    func2taxa["uniref50"] = {}
    func2taxa["ko"] = {}
    func2taxa["module"] = {}
    func2taxa["pathway"] = {}

    # Dict containing level of interest for taxon IDs.
    lineage_ids = {}

    # Line counter.
    lc = 0

    # Read through raw file line by line.
    with open(args.mapping_file[0], "r") as mapping:
        for line in mapping:

            # Skip first line (header).
            if(lc == 0):
                lc += 1
                continue

            # Remove line terminator from end of line.
            line = line.rstrip("\r\n")

            # Split line on whitespace.
            line_split = line.split("\t")

            if line_split[5] not in lineage_ids:
                lineage_ids = {}
                lineage_ids[line_split[5]] = ncbi_taxa_lineage(line_split[5],
                                                               level2keep)

            # Skip because this taxa doesn't have level of interest.
            if(len(lineage_ids[line_split[5]].keys()) == 0):
                continue

            # Add lineage ids to each function in func2taxa dict.
            add_lineage_to_func(d=func2taxa,
                                lin=lineage_ids[line_split[5]],
                                uniref100_id=line_split[7],
                                uniref50_id=line_split[8],
                                uniref90_id=line_split[9],
                                module=line_split[2],
                                pathway=line_split[3],
                                ko=line_split[4],
                                level2keep=level2keep)

    taxa_counts, level_max = get_taxa_counts(d=func2taxa)

    output_num_taxa(counts=taxa_counts, level_max=level_max,
                    level2keep=level2keep)


if __name__ == '__main__':
    main()
