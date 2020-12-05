# TargetedTilingPrimers
Create mutation targeting primers tiling codons of a gene for codon mutagenesis

This repository contains information and a Python script for designing primers to generate mutant libraries of genes with specific desired mutations at each site.

## Codon mutagenesis

Codon mutagenesis can be used to create mutant libraries of a gene with all possible codon mutations.
Such libraries are useful for experiments such as [deep mutational scanning](https://www.ncbi.nlm.nih.gov/pubmed/25075907) or [mutational antigenic profiling](http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006271)

A [script in a separate repository](https://github.com/jbloomlab/CodonTilingPrimers) can be used to generate random mutations at each codon in a gene. This protocol randomly mutates each codon to all other possible codons. This protocol and papers using it are also described in that repository.

Rather than generating mutant libraries with all possible codon mutations at every site in a gene, the script in this repository can be used to generate primers that will result in specific desired mutations at each site in a gene. In this way, a mutant library can be generated which has fewer lethal mutations than random codon mutagenesis.

## Running the script to design primers

The [create_primers.py](create_primers.py) Python script can be used to create primers that tile the codons of a gene in both the forward and reverse direction. Each pair of foward and reverse primers will produce one specified mutation at one specified site. The desired mutations must be given as a csv file argument.

### Required arguments

*seq* : Text file containing the sequence of the gene. The gene itself should be upper case, including the start and stop codons. Any flanking regions should be lower case. There must be >= (primerlength - 3) / 2.0 nucleotides before the first codon mutagenized and after the last codon mutagenized, otherwise there would not be enough flanking nucleotides for a primer to mutate those codons.

*mutations_csv* : csv file containing a table with the mutations the primers will make to seq. The table should have a column, 'site', with site numbers, and a column, 'mutant', with the single letter amino acid mutant to make at that site. Each site can have multiple mutants, or none. The site number should denote the number of the codon in the uppercase gene, with the start codon being 1.

*codon_frequency_csv* : csv file containing a table with codon frequencies used to determine which codons sites will be mutated to. The table should have the columns 'aa', 'codon', and 'frequency'. Each row should have one single letter amino acid in 'aa', one codon corresponding to that amino acid in 'codon', and the frequency of that codon in 'frequency'. The script uses the codon with the highest frequency to replace the codon at a site to make mutations. Codon frequency tables for different organisms can be found [here](https://www.kazusa.or.jp/codon/).

*primerlength* : length of primers. Must be an odd number, so that there is equal length flanking on each side.

*prefix* : string prefix attached to primer names.

The script takes command line arguments; for a listing of how to provide the arguments, type the following to get the help message:

`python create_primers.py -h`

### Optional parameters

There are a variety of optional parameters specifying primer length and melting temperature constraints; the default values for these optional parameters are displayed when you run the program with the `-h` option to get the help message.

## How the script works

The script works as follows:

    1) `mutations_csv` is used to determine what mutations will be made at what sites.

    2) For each site, the frequencies from `codon_frequency_csv` are used to determine what the most frequent codon is for each mutation to be made at that site. This codon will be made by the primer for that mutation at that site.

    3) For each mutant codon, it first makes an ORIGINAL primer of the length specified by `--startprimerlength`

    4) If the original primer has a melting temperature (Tm) greater than the value specified by `--maxprimertm`, then nucleotides are trimmed off one by one (first from the 5' end, then the 3' end, then the 5' end again, etc) until the melting temperature is less than `--maxprimertm` or the length is reduced to `--minlength`.

    5) If the original primer has a Tm greater than `--minprimertm`, then nucleotides are added one-by-one (first to the 3' end, then the 5' end, then the 3' end again, etc) until the melting temperature is greater than `--minprimertm` or the length reaches `--maxlength`.

Note that because the primers are constrained to be between `--minprimerlength` and `--maxprimerlength`, the Tm may not always fall between `--minprimertm` and `--maxprimertm`. This can also happen if a primer initially exceeds `--maxprimertm` but the first trimming that drops it below this value also drops it below `--minprimertm`, or vice-versa if the primer is being extended to increase its melting temperature.

The  *Tm_NN* command of the [MeltingTemp module of Biopython](http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html) is used to calculate Tm of primers.
This calculation is based on nearest neighbor thermodynamics.

The result of running this script is the file specified by `outfile`. It lists the primers.
