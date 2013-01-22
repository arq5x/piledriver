Piledriver.  A no frills, no assumptions pileup engine.


Piledriver is built upon Derek Barnett's excellent BamTools library.
https://github.com/pezmaster31/bamtools


Usage
=====

`piledriver` accepts multiple BAM files at once, and reports the total number
and base quality observed for each allele {A,C,G,T,DEL,INS} among all samples.
In addition, it reports the total number and quality observed for each allele
and sample:

    $ export FTPBASE=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data
    $ bin/bamtools piledriver \
       -in $FTPBASE/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam \
       -in $FTPBASE/HG00111/exome_alignment/HG00111.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam \
    | head
    #chrom	start	end	num_A	num_C	num_G	num_T	num_D	num_I	totQ_A	totQ_C	totQ_G	totQ_T	totQ_D	totQ_I	sample_1	sample_2
    1	9999	10000	1	0	0	0	0	0	2	0	0	0	0	0	0|0,0|0,0|0,0|0	1|2,0|0,0|0,0|0
    1	10000	10001	0	0	0	1	0	0	0	0	0	28	0	0	0|0,0|0,0|0,0|0	0|0,0|0,0|0,1|28
    1	10001	10002	1	0	0	0	0	0	30	0	0	0	0	0	0|0,0|0,0|0,0|0	1|30,0|0,0|0,0|0
    1	10002	10003	2	0	0	0	0	0	32	0	0	0	0	0	1|4,0|0,0|0,0|0	1|28,0|0,0|0,0|0
    1	10003	10004	0	2	0	0	0	0	0	67	0	0	0	0	0|0,1|37,0|0,0|0	0|0,1|30,0|0,0|0
    1	10004	10005	0	2	0	0	0	0	0	66	0	0	0	0	0|0,1|37,0|0,0|0	0|0,1|29,0|0,0|0
    1	10005	10006	0	2	0	0	0	0	0	71	0	0	0	0	0|0,1|38,0|0,0|0	0|0,1|33,0|0,0|0
    1	10006	10007	0	0	0	2	0	0	0	0	0	59	0	0	0|0,0|0,0|0,1|34	0|0,0|0,0|0,1|25
    1	10007	10008	2	0	0	0	0	0	68	0	0	0	0	0	1|39,0|0,0|0,0|0	1|29,0|0,0|0,0|0

The output for each sample is as follows:

`num(A)|totQ(A),num(C)|totQ(C),num(G)|totQ(G),num(T)|totQ(T),num(DEL)|.,num(INS)|totQ(INS)`

For example, 4 `C` alleles were observed for the following sample and the total
quality observed for the 4 alleles was 115.  Deletion alleles have no quality
scores - thus, totQ(DEL) is always set to `.`:

    0|0,4|115,0|0,0|0,0|.,0|0


Installation
============
NOTE: You will need cmake.

1. Clone the repository.
2. `cd piledriver`
3. `mkdir build`
4. `cd build`
5. `cmake ..`
6. `make`


Source
======
All source code for `piledriver` is in `src/toolkit/piledriver*`.  Refinement 
will occur with time.  Suggestions welcome.


To do
=====

The following items must be completed for this to be fully functional.
1. Report stranded counts, etc.
2. Report reference allele.
3. Report uncovered positions.