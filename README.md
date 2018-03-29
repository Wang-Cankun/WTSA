# WTSA


## Usage

This software provides a tool capable of accurately extracting cis-regulatory motifs with their correct lengths even in thousands of promoters from a whole prokaryotic genome. For a fasta file containing promoters, the program outputs a number of profiles, with increasing order of their pvalue, containing the binding sites of specific motifs.

Certain parts of the code uses open-source data structure library codes, including:
- fib <http://resnet.uoregon.edu/~gurney_j/jmpc/fib.html>, copyright information in fib.c
- Mark A. Weiss's data structure codes <http://www.cs.fiu.edu/~weiss/>


## Installation

enter the folder "WTSA" and type "make" then the compiled codes are within the same directory as the source.

## Inputs and outputs

The major program in the provided package is `WTSA`, it can parse standard fasta format of files (allowing squences with different length), and example file are provided. 

To see help and look at all available options.

```console
$ ./WTSA -h (./WTSA)
```

Take a look at `example` (fasta file) first. And try to run WTSA under a specific length, now we can handle motif length more than or equal to 5 (controlled in write_block.c line 256).

```console
$ ./WTSA -i example -l 14
```

For each input file under a specific length l, our program generates a output file, namely, '.closures'file. In '.closures' file, it provides all the closures, representing motif profiles, in the increasing order of their pvalues.

Then try to run WTSA recognizing the correct length in the scope [L,U] by our program automatically

```console
$ ./WTSA -i example -L 14 -U 16
```

L and U are low boundary and up boundary of scope of motif length separately. We use this pair of parameters when we do not know the accurate length in advance. We sort the top n closures under each specific length in the increasing order of their pvalues and save the top o clousres in the ".closures" file. Especially, when the input value of L equals to U, it is equivalent to finding motifs in a specific length. '$ ./WTSA -i example -L 14 -U 14' is equivalent to '$ ./WTSA -i example -l 14'.

## Changelog

### Version 0.1 

Migrated from BoBro2-C-06062017



## Contact

Any questions, problems, bugs are welcome and should be dumped to
Cankun Wang <cankun.wang@sdstate.edu>

