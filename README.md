Piledriver.
===========

**A no frills, no assumptions pileup engine.**


Piledriver is built upon Derek Barnett's excellent BamTools library.
https://github.com/pezmaster31/bamtools


Installation
============
NOTE: You will need cmake.

1. Clone the repository.
2. `cd piledriver`
3. `mkdir build`
4. `cd build`
5. `cmake ..`
6. `make`


Usage
=====

**NOTE**: `piledriver` requires an indexed FASTA file as input.  One can use
`samtools faidx` to index one's FASTA file.

`piledriver` accepts multiple BAM files at once, and reports the total number
and base quality observed for each allele {A,C,G,T,DEL} among all samples.
In addition, it reports the total number and quality observed for each allele
and sample:

    $ export FTPBASE=ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data
	$ export GENOME="hg19.fa"
	$ samtools faidx hg19.fa
    $ bin/bamtools piledriver \
       -fasta hg19.fa
	   -in $FTPBASE/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam \
       -in $FTPBASE/HG00111/exome_alignment/HG00111.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam \
    | head
	#chrom	start	end	ref	depth	r_depth	a_depth	num_A	num_C	num_G	num_T	num_D	num_I	totQ_A	totQ_C	totQ_G	totQ_T	all_ins	sample_1	sample_2
	1	9999	10000	N	1	0	1	1	0	0	0	0	0	2	0	0	0	.	0|0,0|0,0|0,0|0,0|.0|.	1|2,0|0,0|0,0|0,0|.0|.
	1	10000	10001	t	1	1	0	0	0	0	1	0	0	0	0	0	28	.	0|0,0|0,0|0,0|0,0|.0|.	0|0,0|0,0|0,1|28,0|.0|.
	1	10001	10002	a	1	1	0	1	0	0	0	0	0	30	0	0	0	.	0|0,0|0,0|0,0|0,0|.0|.	1|30,0|0,0|0,0|0,0|.0|.
	1	10002	10003	a	2	2	0	2	0	0	0	0	0	32	0	0	0	.	1|4,0|0,0|0,0|0,0|.0|.	1|28,0|0,0|0,0|0,0|.0|.
	1	10003	10004	c	2	2	0	0	2	0	0	0	0	0	67	0	0	.	0|0,1|37,0|0,0|0,0|.0|.	0|0,1|30,0|0,0|0,0|.0|.
	1	10004	10005	c	2	2	0	0	2	0	0	0	0	0	66	0	0	.	0|0,1|37,0|0,0|0,0|.0|.	0|0,1|29,0|0,0|0,0|.0|.
	1	10005	10006	c	2	2	0	0	2	0	0	0	0	0	71	0	0	.	0|0,1|38,0|0,0|0,0|.0|.	0|0,1|33,0|0,0|0,0|.0|.
	1	10006	10007	t	2	2	0	0	0	0	2	0	0	0	0	0	59	.	0|0,0|0,0|0,1|34,0|.0|.	0|0,0|0,0|0,1|25,0|.0|.
	1	10007	10008	a	2	2	0	2	0	0	0	0	0	68	0	0	0	.	1|39,0|0,0|0,0|0,0|.0|.	1|29,0|0,0|0,0|0,0|.0|.

The output for each sample is as follows:

`num(A)|totQ(A),num(C)|totQ(C),num(G)|totQ(G),num(T)|totQ(T),num(DEL)|.,num(INS)|CVS(INS_ALLELES)`

For example, 4 `C` alleles were observed for the following sample and the total
quality observed for the 4 alleles was 115.  Deletion alleles have no quality
scores - thus, totQ(DEL) is always set to `.`:

    0|0,4|115,0|0,0|0,0|.,0|0,3|AA,AA,A


Consider following alignment snapshot:

![](https://raw.github.com/arq5x/piledriver/master/img/igv.png)


Here is the output from underlying BAM file:

    bin/bamtools piledriver -fasta genome.fa -in NA18152.bam | head -40
	#chrom  start   end     ref     depth   r_depth a_depth num_A   num_C   num_G   num_T   num_D   num_I   totQ_A  totQ_C  totQ_G  totQ_T  all_ins sample_1
	chr1	554304	554305	T	5	5	0	0	0	0	5	0	0	0	0	0	117	.	0|0,0|0,0|0,5|117,0|.0|.
	chr1	554305	554306	T	5	5	0	0	0	0	5	0	0	0	0	0	117	.	0|0,0|0,0|0,5|117,0|.0|.
	chr1	554306	554307	T	5	0	5	0	5	0	0	0	0	0	120	0	0	.	0|0,5|120,0|0,0|0,0|.0|.
	chr1	554307	554308	G	5	5	0	0	0	5	0	0	0	0	0	128	0	.	0|0,0|0,5|128,0|0,0|.0|.
	chr1	554308	554309	A	5	5	0	5	0	0	0	0	0	123	0	0	0	.	5|123,0|0,0|0,0|0,0|.0|.
	chr1	554309	554310	C	6	6	0	0	6	0	0	0	0	0	143	0	0	.	0|0,6|143,0|0,0|0,0|.0|.
	chr1	554310	554311	C	6	6	0	0	6	0	0	0	0	0	140	0	0	.	0|0,6|140,0|0,0|0,0|.0|.
	chr1	554311	554312	T	6	6	0	0	0	0	6	0	0	0	0	0	140	.	0|0,0|0,0|0,6|140,0|.0|.
	chr1	554312	554313	T	6	6	0	0	0	0	6	0	0	0	0	0	140	.	0|0,0|0,0|0,6|140,0|.0|.
	chr1	554313	554314	C	6	1	5	0	1	0	0	5	0	0	24	0	0	.	0|0,1|24,0|0,0|0,5|.0|.
	chr1	554314	554315	A	6	0	6	0	0	0	0	6	0	0	0	0	0	.	0|0,0|0,0|0,0|0,6|.0|.
	chr1	554315	554316	G	6	6	0	0	0	6	0	0	0	0	0	143	0	.	0|0,0|0,6|143,0|0,0|.0|.
	chr1	554316	554317	C	6	5	1	0	5	0	0	1	0	0	123	0	0	.	0|0,5|123,0|0,0|0,1|.0|.
	chr1	554317	554318	A	6	1	5	1	0	0	0	5	0	22	0	0	0	.	1|22,0|0,0|0,0|0,5|.0|.
	chr1	554318	554319	A	6	0	6	0	2	0	0	4	0	0	48	0	0	.	0|0,2|48,0|0,0|0,4|.0|.
	chr1	554319	554320	G	6	1	5	0	4	1	1	0	0	0	99	15	22	.	0|0,4|99,1|15,1|22,0|.0|.
	chr1	554320	554321	G	6	6	0	0	0	6	0	0	0	0	0	110	0	.	0|0,0|0,6|110,0|0,0|.0|.
	chr1	554321	554322	T	6	0	6	0	1	0	0	5	0	0	22	0	0	.	0|0,1|22,0|0,0|0,5|.0|.
	chr1	554322	554323	C	6	1	5	0	1	0	0	5	1	0	22	0	0	G	0|0,1|22,0|0,0|0,5|.1|G
	chr1	554323	554324	A	7	7	0	7	0	0	0	0	0	150	0	0	0	.	7|150,0|0,0|0,0|0,0|.0|.
	chr1	554324	554325	A	7	7	0	7	0	0	0	0	0	150	0	0	0	.	7|150,0|0,0|0,0|0,0|.0|.
	chr1	554325	554326	A	7	0	7	0	0	7	0	0	0	0	0	151	0	.	0|0,0|0,7|151,0|0,0|.0|.
	chr1	554326	554327	G	7	7	0	0	0	7	0	0	0	0	0	151	0	.	0|0,0|0,7|151,0|0,0|.0|.
	chr1	554327	554328	G	7	7	0	0	0	7	0	0	0	0	0	153	0	.	0|0,0|0,7|153,0|0,0|.0|.
	chr1	554328	554329	G	7	7	0	0	0	7	0	0	0	0	0	152	0	.	0|0,0|0,7|152,0|0,0|.0|.
	chr1	554329	554330	A	7	7	0	7	0	0	0	0	0	162	0	0	0	.	7|162,0|0,0|0,0|0,0|.0|.
	chr1	554330	554331	G	7	7	0	0	0	7	0	0	1	0	0	156	0	G	0|0,0|0,7|156,0|0,0|.1|G
	chr1	554331	554332	T	7	7	0	0	0	0	7	0	0	0	0	0	158	.	0|0,0|0,0|0,7|158,0|.0|.
	chr1	554332	554333	C	7	7	0	0	7	0	0	0	0	0	156	0	0	.	0|0,7|156,0|0,0|0,0|.0|.
	chr1	554333	554334	C	7	7	0	0	7	0	0	0	0	0	156	0	0	.	0|0,7|156,0|0,0|0,0|.0|.
	chr1	554334	554335	G	8	8	0	0	0	8	0	0	0	0	0	196	0	.	0|0,0|0,8|196,0|0,0|.0|.
	chr1	554335	554336	A	8	8	0	8	0	0	0	0	0	192	0	0	0	.	8|192,0|0,0|0,0|0,0|.0|.
	chr1	554336	554337	A	8	8	0	8	0	0	0	0	1	193	0	0	0	A	8|193,0|0,0|0,0|0,0|.1|A
	chr1	554337	554338	C	9	9	0	0	9	0	0	0	1	0	221	0	0	C	0|0,9|221,0|0,0|0,0|.1|C
	chr1	554338	554339	T	9	9	0	0	0	0	9	0	0	0	0	0	225	.	0|0,0|0,0|0,9|225,0|.0|.
	chr1	554339	554340	A	9	9	0	9	0	0	0	0	0	230	0	0	0	.	9|230,0|0,0|0,0|0,0|.0|.
	chr1	554340	554341	G	9	9	0	0	0	9	0	0	0	0	0	232	0	.	0|0,0|0,9|232,0|0,0|.0|.
	chr1	554341	554342	T	9	9	0	0	0	0	9	0	0	0	0	0	235	.	0|0,0|0,0|0,9|235,0|.0|.
	chr1	554342	554343	C	9	9	0	0	9	0	0	0	0	0	235	0	0	.	0|0,9|235,0|0,0|0,0|.0|.
	chr1	554343	554344	T	9	9	0	0	0	0	9	0	0	0	0	0	232	.	0|0,0|0,0|0,9|232,0|.0|.
	chr1	554344	554345	C	9	9	0	0	9	0	0	0	0	0	230	0	0	.	0|0,9|230,0|0,0|0,0|.0|.
	chr1	554345	554346	A	10	10	0	10	0	0	0	0	0	251	0	0	0	.	10|251,0|0,0|0,0|0,0|.0|.
	chr1	554346	554347	G	10	10	0	0	0	10	0	0	0	0	0	254	0	.	0|0,0|0,10|254,0|0,0|.0|.
	chr1	554347	554348	G	10	10	0	0	0	10	0	0	0	0	0	254	0	.	0|0,0|0,10|254,0|0,0|.0|.
	chr1	554348	554349	C	10	10	0	0	10	0	0	0	0	0	255	0	0	.	0|0,10|255,0|0,0|0,0|.0|.
	chr1	554349	554350	T	10	10	0	0	0	0	10	0	0	0	0	0	252	.	0|0,0|0,0|0,10|252,0|.0|.


Source
======
All source code for `piledriver` is in `src/toolkit/piledriver*`.  The bulk of 
the computation is done in the ``PileDriverPileupFormatVisitor::Visit()``
method. Refinement will occur with time.  Suggestions welcome.


To do
=====

The following items must be completed for this to be fully functional.

1. Report stranded counts, etc.
2. Report uncovered positions.