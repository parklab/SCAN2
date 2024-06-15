BEGIN {
    bulk_first=0;
    sc_first=0;
    fail=0;
    if (length(chrom) == 0) {
        print "chrom not set. use --assign chrom=name";
        fail 1;
        exit 1;
    }
    if (length(max_depth) == 0) {
        print "max_depth not set. use --assign max_depth=N";
        fail=1;
        exit 1;
    }
    if (length(bulk_sample) == 0) {
        print "bulk_sample not set. use --assign bulk_sample=name";
        fail=1;
        exit 1;
    }
    if (length(sc_sample) == 0) {
        print "sc_sample not set. use --assign sc_sample=name";
        fail=1;
        exit 1;
    }
    for (i = 0; i <= max_depth; ++i) {
        for (j = 0; j <= max_depth; ++j) {
            table[i][j] = 0;
        }
    }
}
NR==1 {
    # Ensure the named bulk and single cell samples are present
    if ($1 != bulk_sample && $2 != bulk_sample) {
        print "bulk sample " bulk_sample " is not present in the input file header, exiting";
        fail=1;
        exit 1;
    } else {
        if ($1 == bulk_sample) {
            bulk_first = 1;
        }
    }
            
    if ($1 != sc_sample && $2 != sc_sample) {
        print "sc sample " sc_sample " is not present in the input file header, exiting";
        fail=1;
        exit 1;
    } else {
        if ($1 == sc_sample) {
            sc_first = 1;
        }
    }
}
NR>1 {
    count1 = $1 <= max_depth ? $1 : max_depth;
    count2 = $2 <= max_depth ? $2 : max_depth;
    table[count1][count2] = table[count1][count2] + 1;
}
END {
    if (fail == 0) {
        print "chrom=" chrom "\tsc_sample=" sc_sample "\tbulk_sample="bulk_sample "\tsc_max_depth=" max_depth "\tbulk_max_depth=" max_depth;
        # Expected format is rows = depth of single cell, columns = depth of bulk
        if (sc_first == 1) {
            for (i = 0; i <= max_depth; ++i) {
                printf "%ld", table[i][0];
                for (j = 1; j <= max_depth; ++j) {
                    printf "\t%ld", table[i][j];
                }
                printf "\n";
            }
        } else {
            for (i = 0; i <= max_depth; ++i) {
                printf "%ld", table[0][i];
                for (j = 1; j <= max_depth; ++j) {
                    printf "\t%ld", table[j][i];
                }
                printf "\n";
            }
        }
    }
}
