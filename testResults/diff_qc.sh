#!/usr/bin/env bash
#
# quick diff of NEW vs OLD
#
echo "< NEW QC.regression.new.tsv"
echo "> OLD QC.regression.tsv"
diff -w <(cut -f 5-7,9-10 QC.regression.new.tsv) <(cut -f 5-7,9-10 QC.regression.tsv)
