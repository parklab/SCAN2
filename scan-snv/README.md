Inputs to SCAN-SNV are:
  * Two GATK tables produced by the totab awk script. One table should be
    generated with stringent read mapping quality requirements (e.g., with
    -mmq 60) and the other should be very permissive (e.g., -mmq 1).
  * R data files containing objects:
    * data          - hSNP training data. Locations, phases and numbers of
                      alt and ref reads in the single cell being analyzed.
    * fits          - covariance function parameter fits from gridfit
  * A desired false discovery rate. There is no guarantee that this rate will
    be met.

SCAN-SNV creates intermediary files containing:
  * somatic.ab      - AB predictions at all somatic SNV candidate loci.
                      Predictions are in the form of parameters to the
                      normal posterior distribution.
  * hsnp.cigars     - A table in which each row represents one hSNP and gives
                      counts of several CIGAR operations in reads covering the
                      hSNP.
  * somatic.cigars  - Same as above, but for somatic candidate sites.


SCAN-SNV output:
  * somatic_gt.rda  - R data file containing an object named "gt".
                      gt$somatic[gt$somatic$pass,] is the set of called sSNVs.
