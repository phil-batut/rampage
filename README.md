# rampage repository
Analysis scripts for RAMPAGE data (Batut et al., Genome Res. 2013) 



# rampagePeakCaller.py
# Philippe BATUT
# Gingeras Lab, CSHL
# Sept. 22, 2011





This script takes as input a WIGGLE file describing the density of RAMPAGE read 5' ends
over the genome. For a sliding window of size k sliding through each position across the
genome, it computes the probability of obtaining the observed number of RAMPAGE tags under
the null model. The null model is a negative binomial distribution whose mean is the
observed average density of tags per window, and whose dispersion parameter is user-defined.

ptionally, it can also be provided with a second wiggle file describing the density of
RAMPAGE downstream reads (i.e., the read of the mate pair that starts at the RT priming
site). The downstream read coverage is used to subtract a pseudocount attributable to
background signal prior to computing the p-value for each window. The weight to be given
to this backgroung estimate is a user-defined parameter.

The p-values are corrected by the method of Benjamini & Hochberg and thresholded to allow
for an FDR below a value defined by the user. All significant windows within a distance d
of each other are then merged, which generates the final peaks.


INPUT:       - Argument 1:            File providing the size of each chromosome in the reference
                                      genome (Format <chrName\tchrSize\n>).
             - Argument 2:            RAMPAGE tag density file (Wiggle).
                                      The track line for each strand is required to contain a name
                                      field ending with "_+" or "_-", from which the strand is
                                      deduced.


OPTIONS:     --width (-w):            Sliding window width. Default: 15 bases.
             --NBinom (-n):           Dispersion parameter of the negative binomial. Default: 0.6.
             --read2-background:      Downstream read density file (wiggle).
             --bg-weight (-w)         Weight to be given to the background signal. Default: 0.5.
             --fdr (-f):              False discovery rate. Default: 10^(-6).
             --exclude (-x):          Coma-separated list of chromosomes to exclude. Default: "chrM".
             --distance (-d):         Maximum distance for merging neighboring windows into a single peak.
             --threads (-t):          Number of threads to be used. Default: 1.


In our experience, the dispersion parameter (--NBinom) and the FDR (--fdr) are the most critical
ones to change when optimizing the specificity of the peakcalling.  
