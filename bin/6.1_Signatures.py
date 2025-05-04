#!/usr/bin/env python3

#input (maf, vcf..)

#Format: 
#Hugo_Symbol	Chromosome	Start_Position	Variant_Type	Reference_AllelTumor_Seq_Allele2	Tumor_Sample_Barcode
#Gene1	1	158862794	SNP	G	A	P01
#Gene2	19	58864867	SNP	A	G	P02


from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
import sigProfilerPlotting as sigPlt
from SigProfilerAssignment import Analyzer as Analyze


matrices = matGen.SigProfilerMatrixGeneratorFunc("luad", "GRCh37", \
                                                 "/Signatures/", \
                                                    plot=True, exome=True)

sigPlt.plotSBS("/Signatures/output/SBS/luad.SBS96.exome", \
                "/Signatures/output/SBS/", "luad", "96", \
                    percentage=False)

Analyze.cosmic_fit("/Signatures/output/SBS/luad.SBS96.exome", \
                   "/labs/lgfc/signatures/SPA/", \
                    input_type="matrix", \
                    context_type="96", \
                    cosmic_version=3.4, \
                    exome=True, \
                    genome_build="GRCh37", \
                    make_plots=True)
