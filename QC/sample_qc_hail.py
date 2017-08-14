import numpy as np
import pandas as pd
from collections import Counter
from math import log, isnan
from pprint import pprint
from hail import *

hc = HailContext()
vds = hc.read("gs://gnomad-berylc/MYOSEQ.vds/")
vds = vds.split_multi()

table = hc.import_table('gs://gnomad-berylc/MYOSEQ.annotations.tsv', impute= True).key_by('Sample')
vds = vds.annotate_samples_table(table, root='sa')

#Sample and variant QC first
vds = vds.sample_qc()
vds = vds.variant_qc().cache()

#Filter to biallelic SNPs with AF > 0.01
vds_sex_check = vds.filter_multi()
vds_sex_check = vds_sex_check.filter_variants_expr('va.filters.isEmpty() && va.qc.AF > 0.01 && v.altAllele.isSNP()', keep=True)

#Impute sex 
vds_sex_check = vds_sex_check.impute_sex()

#Also mimic gnomad sex check 
#X heterozygosity
vds = vds.annotate_samples_expr('sa.NVarsChrX = gs.filter(g => v.contig == "X").count()')
vds = vds.annotate_samples_expr('sa.NHetVarsChrX = gs.filter(g => v.contig == "X" && g.isHet()).count()')

#Chrom 20 coverage
vds = vds.annotate_samples_expr('sa.Cov20 = gs.filter(g => v.contig == "20").map(g => g.dp).stats().sum')

#Chrom Y coverage
vds_y = vds.filter_variants_expr('v.contig == "Y"', keep = True)
intervals = map(Interval.parse, ['Y:10001-2649520', 'Y:59034050-59363566'])
vds_y = vds_y.filter_intervals(intervals, keep = False)
vds_y = vds_y.annotate_samples_expr('sa.CovY = gs.filter(g => v.contig == "Y").map(g => g.dp).stats().sum')

vds_y.count_variants()

#Write file
vds.export_samples("gs://gnomad-berylc/output/sample_qc_081417.tsv", 'Sample = s, sa.qc.*, sa.imputesex.*,  CovY =  sa.CovY, Cov20 = sa.Cov20, NHetX = sa.NHetVarsChrX, NTotalX = sa.NVarsChrX')

#Export metrics to plot
#vds.export_samples("gs://gnomad-berylc/output/sample_qc.tsv", 'Sample = s, sa.qc.*, sa.imputesex.*')

#Filter based on QC metrics
vds = vds.filter_samples_expr('sa.qc.dpMean > 15 && sa.qc.callRate > 0.9 && sa.qc.gqMean > 50 && sa.qc.rTiTv < 2.65 && sa.qc.rHetHomVar >1 && sa.qc.rInsertionDeletion < 1.2')
#This leaves 14,683 samples

#Filter based on sex imputation
vds = vds.filter_samples_expr('sa.imputesex.Fstat > -0.3', keep = True)
vds = vds.filter_samples_expr('sa.imputesex.Fstat > 0.35 && sa.imputesex.Fstat < 0.5', keep = False)
#First filters 106 samples, second filters 6, can make these stricter in later iterations

#Filter three additional MyoSeq samples where the sex in ped does not match inferred sex
sex_in_ped_inconsistent = ["MKHA009", "MTEH052","MVAL013","MVAL012"]
vds = vds.filter_samples_list(sex_in_ped_inconsistent, keep = False)

#This leaves 14,567 samples

#Calculate relatedness on autosomal biallelic variants w/ > 99% call rate > 1% AF. Since we already split multiallelic, we need to remove split variants
#vds_gnomad_filters = vds.filter_variants_expr('va.wasSplit', keep = False)
#vds_gnomad_filters = vds_gnomad_filters.filter_variants_expr('v.contig != "X" && v.contig != "Y" && v.contig != "MT" ')
#vds_gnomad_filters = vds.variant_qc().cache()

#vds_gnomad_filters = vds_gnomad_filters.filter_variants_expr('va.qc.callRate > 0.99 && va.qc.AF > 0.01' , keep = True).ld_prune(r2 = 0.1)
#ibd_calculations = vds_gnomad_filters.ibd(min = 0.2)
#ibd_calculations_gnomad.export("gs://gnomad-berylc/output/ibd.calculations.on.filtered.data.gnomadfilters.072517.tsv”)

#Removing samples manually
to_remove_relatedness = hc.import_table("gs://gnomad-berylc/samples_filter_relatedness.after.initial.qc.072517.txt", no_header=True).key_by('f0')
vds = vds.filter_samples_table(to_remove_relatedness, keep = False)
#This leaves 14,285 samples
#vds.write("gs://gnomad-berylc/filtered.MYOSEQ.sex.qc.relatedness.vds")

#Filter to only European individuals 
to_keep_europeans = hc.import_table("gs://gnomad-berylc/samples_keep_european.tsv", no_header=True).key_by("f0")
vds = vds.filter_samples_table(to_keep_europeans, keep = True)
#This leaves 13,909 samples
