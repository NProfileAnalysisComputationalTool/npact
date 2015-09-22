# About #

NPACT (N-Profile Analysis and Computation Tool) implements methods to
identify in nucleotide sequences regions of statistically-significant
compositional three-base periodicity associated to ORF structures,
to compare the newly-identified ORFs with pre-annotated collections
of genes and to highlight potential new genes not included in the annotation.
The primary output is a graphical representation of sequence features
represented by frame-specific compositional profiles (by default, GC content)
and sequence segments with significant periodicity (Hits) together with
representation of the pre-annotated genes and, in a separate track, newly
identified ORFs not included in the annotation.

NPACT accepts input sequences of any length and its dynamical graphical
representations allow easy browsing and analysis of an entire prokaryotic genome.

This NPACT project was developed in 2011-2015 by Acceleration.net under
contract from Luciano Brocchieri at the [University of Florida
Genetics Institute](http://www.ufgi.ufl.edu/) based on sequence-analysis
methods and computational tools developed by Luciano Brocchieri at the
[University of Florida Genetics Institute](http://www.ufgi.ufl.edu/).


# The Java branch #

The code in here was an attempt to make a java frontend for NPACT.

I believe it mostly worked but has gotten dated and at this point is
included only for posterity.

# Structure #

To be built it needed both the compiled C programs and the java jar of
these files. The directory structure of the linux build is as follows:

- NPACT.jar
- analysis/
- tables/

## Analysis folder ##

The analysis folder holds the compiled C programs (`nprofile`,
`Allplots`, etc). We currrently have these built at
`pynpact/pynpact/bin/`.

## Tables folder ##

This holds the scoring tables used by `acgt_gamma`. We currently store
these at `pynpact/pynpact/data`
