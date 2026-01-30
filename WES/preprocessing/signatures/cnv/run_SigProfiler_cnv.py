from SigProfilerAssignment import Analyzer as Analyze

Analyze.cosmic_fit('/michorlab/jacobg/Ellisen/signatures/SigProfiler/all_cnvs.tsv', 
                   '/michorlab/jacobg/Ellisen/signatures/SigProfiler/cnv_results', input_type="seg:ASCAT_NGS", 
                   context_type="96", collapse_to_SBS96=False, cosmic_version=3.4, 
                   exome=True, genome_build="GRCh38", signature_database=None,
                   export_probabilities=True,
                   export_probabilities_per_mutation=True, make_plots=True,
                   sample_reconstruction_plots='both', verbose=True)


