"""Script for creating gene assembly primers.

Primers containing specified mutageneic codons at their center are tiled along a gene.

Edited by Caelan Radford in 2020, from the CodonTilingPrimers code written by Jesse Bloom in 2013 and edited by Adam Dingens Nov 2015 to generate primers of differing lengths to all have a Tm of ~60C.

This script first makes an ORIGINAL primer of specified length (default 37 bps).
If the ORIGINAL primer has a Tm of greater than MaxPrimerTm, then nucleotides are trimmed off one by one (first 5', then 3', then 5' etc) until the Tm is less than MaxPrimerTm. Note that this could be over a degree less than the MaxPrimerTm.
If the ORIGINAL primer has a Tm less than MinPrimerTm, then nucelotides are added one by one (first 3', then 5', then 3' etc) until the Tm is over MinPrimerTm. Note that this could be over a degree more than the MinPrimerTm
If the ORIGINAL primer has a Tm of less than MaxPrimerTm but greater than MinPrimerTm, it is not altered.
The primers are constrained to be between MinPrimerlength and MaxPrimerLength bps long. The Tm of some MaxPrimerLength primers may not be > MinPrimerTemp, and the Tm of some MinPrimerLength primers may not be < MaxPrimerTm.

For command line arguments, run::

    python create_primers.py -h

The  Tm_NN command of the MeltingTemp Module of Biopython (http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html) is used to calculate Tm of primers.
This calculation is based on nearest neighbor thermodynamics. nucelotides labeled N are given average values in the Tm calculation.
It is possible to vary salt concentration and other addatives if needed."""


import os
import sys
import math
import re
import argparse
import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq


def Parser():
    """Returns command line parser."""
    parser = argparse.ArgumentParser(
            description='Script by Caelan Radford, Adam Dingens, and Jesse Bloom to design specified codon tiling primers with specific melting temperatures.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            )

    parser.add_argument('sequencefile', help="the name of a file giving the sequence for which we are designing the primers. This file should only contain the sequence, and should not have any headers or other content. For the sequence, make the 5' and 3' ends that you do not want to mutate in lower case. Make the portion of the coding sequence that you want to tile with ambiguous codons in upper case. Typically, for example, the first site you mutate would be the initial start codon, so the first *ATG* would be the first upper case letter. You must have at least *(startprimerlength - 3) / 2* nucleotides in lower case at each end of the upper case sequence that you are mutating. This is because at least this much flanking sequence is needed to design primers of the indicated length; more sequence may be required if the primer at either end is extended beyond the startprimerlength.")
    parser.add_argument('mutations_csv', help="the name of a csv file containing a table of specific mutations per site you want to make. It should have one column ('site') with site numbers and one column ('mutant') with amino acid mutants to make at that site. Site 1 is the first codon in the uppercase sequence.")
    parser.add_argument('codon_frequency_csv', help="the name of a csv file with a table of codon frequencies to be used to determine which codons to mutate to. It should have one column 'aa' with single letter amino acids, one column 'codon' with codons, and one column 'frequency' with codon frequencies. The highest frequency codon for each amino acid is used.")
    parser.add_argument('primerprefix', help="prefix name to be added to each primer")
    parser.add_argument('outfile', help='name of primer output file')
    parser.add_argument('--startprimerlength', type=int, help='starting primer length', default=37)
    parser.add_argument('--maxprimertm', type=float, help="Upper temperature limit for primers.", default=61)
    parser.add_argument('--minprimertm', type=float, help="Lower temperature limit for primers.", default=60)
    parser.add_argument('--minlength', type=int, help='Minimum primer length', default=25)
    parser.add_argument('--maxlength', type=int, help='Maximum primer length', default=51)
    return parser


def ReverseComplement(seq):
    """Returns reverse complement of sequence. Preserves upper/lower case."""
    d = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N', 'n':'n', 'S':'S', 's':'s', 'K':'M', 'k':'m'}
    rc = [d[nt] for nt in seq]
    rc.reverse()
    return ''.join(rc)


def CreateMutForOligosVarLength(seq, mutations_csv, primerlength, prefix, maxprimertm, minprimertm, maxlength, minlength):
    """Creates oligos to tile a gene and introduce specified amino acid mutations at each site.

    *seq* : sequence of the gene. The gene itself should be upper case,
    including the start and stop codons. Any flanking regions should be lower
    case.
    There must be >= (primerlength - 3) / 2.0 nucleotides before the first codon
    mutagenized and after the last codon mutagenized.

    *mutations_csv* : csv file containing a table with the mutations the primers
    will make to seq. The table should have one column, 'site', with site
    numbers, and one column, 'mutant', with the single letter amino acid mutant
    to make at that site. Each site can have mulitple mutants, or none. The site
    number should denote the number of the codon in the uppercase gene, with the
    start codon being 1.

    *primerlength* : length of primers. Must be an odd number, so that equal length
    flanking on each side.

    *prefix* : string prefix attached to primer names.

    Tiles primers across the gene in the forward direction. A unique primer is
    made for each mutation specified by mutations_csv at each site.
    Primers are named as follows:

    f"{prefix}-for-mut{i}{mutant}" for 5' tiling primers, where i is the site
    mutagenized in the gene and mutant is the amino acid the site is mutated to
    using that primer.

    Returns a list of all these primers as *(name, sequence)* 2-tuples.
    """
    n = len(seq)
    assert primerlength % 2 == 1, "primer length not odd"
    initial_flanklength = (primerlength - 3) // 2
    upperseq = ''.join([nt for nt in seq if nt.istitle()])
    assert upperseq in seq, "upper case nucleotides not substring"
    assert len(upperseq) % 3 == 0, "length of upper case not multiple of 3"
    ncodons = len(upperseq) // 3

    # Read in the mutations csv and make sure it is the right format
    df = pd.read_csv(mutations_csv)
    try:
        sites = df['site'].tolist()
    except:
        print('mutations_csv must have column "site"')
    try:
        mutations = df['mutant'].tolist()
    except:
        print('mutations_csv must have column "mutant"')
    mutations_dict = {}
    for i, site in enumerate(sites):
        if site not in mutations_dict:
            mutations_dict[site] = [mutations[i]]
        else:
            mutations_dict[site].append(mutations[i])
    # Make sure this is enough flanking sequence to make primers
    initial_flanklength = (primerlength - 3) // 2
    firstcodon = min(sites)
    lastcodon = max(sites)
    startupper = seq.index(upperseq) + ((firstcodon - 1) * 3)
    if startupper < initial_flanklength:
        raise ValueError("not enough 5' flanking nucleotides")
    if n - len(upperseq) - seq.index(upperseq) + (len(upperseq) - (lastcodon - 1) * 3) < initial_flanklength:
        raise ValueError("not enough 3' flanking nucleotides")
    # Read in the codon frequency table, make dict for back translating
    df = pd.read_csv('codon_frequencies.csv')
    aas = df['aa'].tolist()
    codons = df['codon'].tolist()
    frequencies = df['frequency'].tolist()
    back_t_dict = {}
    for i, codon in enumerate(codons):
        if aas[i] not in back_t_dict:
            back_t_dict[aas[i]] = {}
            back_t_dict[aas[i]][frequencies[i]] = codon
        else:
            back_t_dict[aas[i]][frequencies[i]] = codon
    # Iterate through the codons and make the specific primers, with most
    # frequent codon
    primers = []
    for icodon in range(ncodons):
        # Determine what the insert should be for the codon
        # Replace ambiguous_codon with codon insert later
        codon_inserts = []
        if icodon + 1 in mutations_dict:
            for mutation in mutations_dict[icodon + 1]:
                max_frequency = max(back_t_dict[mutation])
                codon_insert = back_t_dict[mutation][max_frequency]
                i = startupper + icodon * 3
                primer = "%s%s%s" % (seq[i - initial_flanklength : i], codon_insert, seq[i + 3 : i + 3 + initial_flanklength])
                name = f"{prefix}-for-mut{icodon + 1}{mutation}"
                primerseq = Seq(primer)

                TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False))
                add_3 = True
                minus_5 = True
                flank5 = flank3 = initial_flanklength
                if float(TmNN) > float(maxprimertm):
                    while float(TmNN) > float(maxprimertm) and len(primer) > minlength:
                        if minus_5:
                            flank5 -= 1
                            primer = "%s%s%s" % (seq[i - (flank5) : i], codon_insert, seq[i + 3 : i + 3 + flank3])
                            minus_5 = False
                        else:
                            flank3 -= 1
                            primer =  "%s%s%s" % (seq[i - (flank5) : i], codon_insert, seq[i + 3 : i + 3 + flank3])

                            minus_5 = True
                        if seq.index(upperseq) + ((firstcodon - 1) * 3) < flank5:
                            raise ValueError("not enough 5' lower case flanking nucleotides")
                        if n - len(upperseq) - seq.index(upperseq) + (len(upperseq) - (lastcodon - 1) * 3) < flank3:
                            raise ValueError("not enough 3' lower case flanking nucleotides")


                        primerseq = Seq(primer)
                        TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False))
                        primerlength= len(primer)
                else:
                    if float(TmNN) < float(minprimertm):
                        while float(TmNN) < float(minprimertm) and len(primer) < maxlength:
                            if add_3:
                                flank3 += 1
                                primer = "%s%s%s" % (seq[i - (flank5) : i], codon_insert, seq[i + 3 : i + 3 + flank3])
                                add_3 = False
                            else:
                                flank5 +=1
                                primer = "%s%s%s" % (seq[i - (flank5) : i], codon_insert, seq[i + 3 : i + 3 + flank3])
                                add_3 = True
                            primerseq = Seq(primer)
                            TmNN = ('%0.2f' % mt.Tm_NN(primerseq, strict=False))
                            primerlength= len(primer)
                            if seq.index(upperseq) + ((firstcodon - 1) * 3) < flank5:
                                raise ValueError("not enough 5' lower case flanking nucleotides")
                            if n - len(upperseq) - seq.index(upperseq) + (len(upperseq) - (lastcodon - 1) * 3) < flank3:
                                raise ValueError("not enough 3' lower case flanking nucleotides")

                    else:
                        pass
                primers.append((name, primer))
    #print primers
    return primers


def main():
    parser = Parser()
    args = vars(parser.parse_args())

    print("Read the following command line arguments")
    for (argname, argvalue) in args.items():
        print("\t%s = %s" % (argname, argvalue))


    primerlength = args['startprimerlength']

    if (primerlength <=3 ) or (primerlength % 2 == 0):
        raise ValueError("Does not appear to be valid primer length: %d" % primerlength)

    sequencefile = args['sequencefile']
    if not os.path.isfile(sequencefile):
        raise IOError("Cannot find sequencefile %s" % sequencefile)
    sequence = open(sequencefile).read()
    sequence = sequence.replace(' ', '')
    sequence = sequence.replace('\n', '')
    print("Read a sequence of length %d from %s:\n%s" % (len(sequence), sequencefile, sequence))
    outfile = args['outfile']
    primerprefix = args['primerprefix']
    print("The primers will be named with the prefix %s" % (primerprefix))

    mutations_csv = args['mutations_csv']
    if not os.path.isfile(mutations_csv):
        raise IOError(f"Cannot find mutations_csv {mutations_csv}")

    # Design forward mutation primers
    mutforprimers = CreateMutForOligosVarLength(sequence, mutations_csv, primerlength, primerprefix, args['maxprimertm'], args['minprimertm'], args['maxlength'], args['minlength'])
    print("Designed %d mutation forward primers." % len(mutforprimers))

    # Design reverse mutation primers
    mutrevprimers = [(name.replace('for', 'rev'), ReverseComplement(seq)) for (name, seq) in mutforprimers]
    print("Designed %d mutation reverse primers." % len(mutrevprimers))

    # Print out all of the primers
    primers = mutforprimers + mutrevprimers
    print("This gives a total of %d primers." % len(primers))

    print("\nNow writing these primers to %s" % outfile)
    iplate = 1
    f = open(outfile, 'w')
    for primers in [mutforprimers, mutrevprimers]:
        f.write("\r\nPlate %d\r\n" % iplate)
        n_in_plate = 0
        for (name, primer) in primers:
            f.write("%s, %s\r\n" % (name, primer))
            n_in_plate += 1
            if n_in_plate == 96:
                n_in_plate = 0
                iplate += 1
                f.write("\r\nPlate %d\r\n" % iplate)
        if n_in_plate:
            iplate += 1



main() # run the main program
