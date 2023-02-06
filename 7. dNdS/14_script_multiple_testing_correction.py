#!/usr/bin/env python

import sys
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
import numpy as np

dnds_out = pd.read_csv(sys.argv[1], sep = '\t')
dnds_out = pd.DataFrame(dnds_out)
dnds_out["pvalue"] = pd.to_numeric(dnds_out["pvalue"])
FDR_corr = float(sys.argv[2])

outputfilename=sys.argv[1].replace(".output",f'.multipletest-correction.FDR{FDR_corr}.output')
output=open(outputfilename,'w')

# BH correction
rejected, p_adjusted = fdrcorrection(dnds_out["pvalue"], alpha=FDR_corr, is_sorted=False)
dnds_out[f'FDR{FDR_corr}'] = rejected.tolist()

# print stats on multiple testing
sign = dnds_out[f'FDR{FDR_corr}'][dnds_out[f'FDR{FDR_corr}'] == True].count()
total = len(dnds_out)
perc = round(int(sign)/int(total)*100, 2)
print(f"{sign} of {total} tests ({perc} %) passed multiple-testing correction (BH)")

dnds_out.to_csv(outputfilename, index = False, header = True, sep = '\t', mode = 'w', na_rep = 'NA', float_format='{:f}'.format)
