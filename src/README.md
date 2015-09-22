# Info
The code in here was an attempt to make a java frontend for NPACT.

I believe it mostly worked but has gotten dated and at this point is
included only for posterity.

# Structure

To be built it needed both the compiled C programs and the java jar of
these files. The directory structure of the linux build is as follows:

- NPACT.jar
- analysis/
- tables/

## Analysis folder

The analysis folder holds the compiled C programs (`nprofile`,
`Allplots`, etc). We currrently have these built at
`pynpact/pynpact/bin/`.

## Tables folder

This holds the scoring tables used by `acgt_gamma`. We currently store
these at `pynpact/pynpact/data`
