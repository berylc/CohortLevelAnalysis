import numpy as np
import pandas as pd
from collections import Counter
from math import log, isnan
from pprint import pprint
from hail import *


def read_in_vds(f):
	vds = hc.read(f)
	vds = vds.split_multi()
	return(vds)


def filter_vds_variants(v, ibd=False):
	v = v.filter_variants_expr('va.wasSplit', keep=False)
	v = v.filter_variants_expr(
		'va.filters.isEmpty() &&'
		'va.qc.AF > 0.01 && v.altAllele().isSNP()', keep=True)

	if ibd:
		v = v.filter_variants_expr(
			'v.contig != "X" && v.contig != "Y" &&'
			'v.contig != "MT" ')

		v = v.filter_variants_expr(
			'va.qc.callRate > 0.99', keep=True).ld_prune(r2=0.1)

	return(v)


def filter_samples_qc_metrics(
	v, max_mean_dp=120, min_mean_dp=20, min_callrate=0.95,
	min_mean_gq=60, max_titv=2.6, min_het_homvar=1.2, max_ins_del=1):

	v = v.filter_samples_expr(
		'sa.qc.dpMean > %(min_mean_dp)s &&'
		'sa.qc.dpMean < %(max_mean_dp)s &&'
		'sa.qc.callRate > %(min_callrate)s &&'
		'sa.qc.gqMean > %(min_mean_gq)s &&'
		'sa.qc.rTiTv < %(max_titv)s &&'
		'sa.qc.rHetHomVar > %(min_het_homvar)s &&'
		'sa.qc.rInsertionDeletion < %(max_ins_del)s' % locals())

	return(v)


def get_samples_to_filter_sex(
	v, Fstat_male=0.9, Fstat_female_min=-0.25, Fstat_female_max=0.3):

	v1 = v.filter_samples_expr(
		'sa.imputesex.Fstat > %(Fstat_female_max)s &&'
		'sa.imputesex.Fstat < %(Fstat_male)s ' % locals(), keep=True).sample_ids

	v2 = v.filter_samples_expr(
		'sa.imputesex.Fstat < %(Fstat_female_min)s' % locals(), keep=True).sample_ids

	samples_to_filter = v1 + v2 + ["MKHA004"]

	return(samples_to_filter)


def add_klinfelter_check_information_and_write(v, date):
	v = v.annotate_samples_expr(
		'sa.NVarsChrX = gs.filter(g => v.contig == "X").count(),'
		'sa.NHetVarsChrX = gs.filter(g => v.contig == "X" &&'
		'g.isHet()).count(), sa.Cov20 = '
		'gs.filter(g => v.contig == "20").map(g => g.dp).stats().sum')

	vds_y = v.filter_variants_expr('v.contig == "Y"', keep=True)
	intervals = map(Interval.parse, ['Y:10001-2649520', 'Y:59034050-59363566'])
	vds_y = vds_y.filter_intervals(intervals, keep=False)

	vds_y.export_samples(
		"gs://gnomad-berylc/output/sample_qc_%s.tsv",
		'Sample = s, sa.qc.*, sa.imputesex.*,'
		'CovY =  sa.CovY, Cov20 = sa.Cov20, NHetX = sa.NHetVarsChrX,'
		'NTotalX = sa.NVarsChrX' % date)


def write_relatedness_calculations(v, date):
	vds_ibd = filter_vds_variants(v, ibd=True)
	ibd_calculations = vds_ibd.ibd(min=0.2)
	ibd_calculations.export("gs://gnomad-berylc/output/ibd.calculations.on.filtered.data.gnomadfilters.%s.tsv" % (date))


def write_PCA(v, date, only_known_pop=False):
	gnomad_pca_vds = hc.read("gs://gnomad-genomes/sampleqc/gnomad.pca.vds/")

	# Remove samples in gnomad PCA that overlap with the myoseq callset
	mapping_dict = {sample: sample.replace("exome_", "") for sample in gnomad_pca_vds.sample_ids}
	gnomad_pca_vds = gnomad_pca_vds.rename_samples(mapping_dict)
	gnomad_pca_vds = gnomad_pca_vds.filter_samples_list(vds.sample_ids, keep=False)

	if only_known_pop:
		gnomad_known_pop_samples = hc.import_table(
			"gs://gnomad-berylc/gnomad.known.population.sample.ids.list",
			no_header=True).key_by('f0')

		gnomad_pca_vds = gnomad_pca_vds.filter_samples_table(gnomad_known_pop_samples, keep=True)
		date = "only.known.gnomad.pop.myoseq.%s" % date

	# Harmonize samples schemas of the two vds
	gnomad_pca_vds = gnomad_pca_vds.annotate_samples_expr('sa={}')
	v = v.annotate_samples_expr('sa={}')

	combined_vds = v.join(gnomad_pca_vds)
	combined_vds = combined_vds.variant_qc().cache()  # 146966 samples
	combined_vds = filter_vds_variants(combined_vds, ibd=True)

	pca = combined_vds.pca('sa.pca', k=7)

	pca.export_samples("gs://gnomad-berylc/output/filtered.MYOSEQ.%s.tsv" % (date), 'Sample=s, sa.pca.*')


def myoseqPCA(v, date):
	v = v.filter_variants_expr('va.wasSplit', keep=False)
	v = v.filter_variants_expr(
		'va.filters.isEmpty() && v.altAllele().isSNP()', keep=True)
	v = v.filter_variants_expr(
		'v.contig != "X" && v.contig != "Y" && v.contig != "MT" ')
	v = v.filter_samples_expr("sa.IsMyoseq").variant_qc()
	v = vds_for_subpop_first_filter.filter_variants_expr(
		'va.qc.callRate > 0.99 && va.qc.AF > 0.001', keep=True)
	v01 = v.ld_prune(r2=0.1)
	v05 = v.ld_prune(r2=0.5)
	v09 = v.ld_prune(r2=0.9)

	subpop_pca = vds_subpop.pca('sa.pca', k=7)


def main():
	hc = HailContext()

	vds = read_in_vds("gs://gnomad-berylc/MYOSEQ.vds/")

	table = hc.import_table(
		'gs://gnomad-berylc/MYOSEQ.annotations.tsv',
		impute=True).key_by('Sample')
	vds = vds.annotate_samples_table(table, root='sa')

	vds = vds.sample_qc()
	vds = vds.variant_qc().cache()

	vds = filter_samples_qc_metrics(vds)

	# Sex check
	sex_check_vds = filter_vds_variants(vds)
	sex_check_vds = sex_check_vds.impute_sex()
	sex_check_samples_to_filter = get_samples_to_filter_sex(sex_check_vds)
	# add_klinfelter_check_information_and_write(sex_check_vds, "081417")

	vds = vds.filter_samples_list(sex_check_samples_to_filter, keep=False)

	# Relatedness check
	# write_relatedness_calculations(vds, "081417")

	to_remove_relatedness = hc.import_table("gs://gnomad-berylc/samples_filter_relatedness.after.initial.qc.081417.txt", no_header=True).key_by('f0')
	vds = vds.filter_samples_table(to_remove_relatedness, keep=False)

	# PCA
	# write_PCA(vds, "082517")
	to_keep_population = hc.import_table("gs://gnomad-berylc/output/samples_keep_pca.after.qc.081517.txt", no_header=True).key_by('f0')
	vds = vds.filter_samples_table(to_keep_population, keep=True)


if __name__ == "__main__":
	hc = HailContext()
