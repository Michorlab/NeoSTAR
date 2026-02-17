from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit('/michorlab/jacobg/Ellisen/signatures/SigProfiler/input_mafs', 
                   '/michorlab/jacobg/Ellisen/signatures/SigProfiler/results', input_type="vcf", 
                   context_type="96", collapse_to_SBS96=True, cosmic_version=3.4, 
                   exome=True, genome_build="GRCh38", signature_database=None,
                   exclude_signature_subgroups=['Artifact_signatures'], export_probabilities=True,
                   export_probabilities_per_mutation=True, make_plots=True,
                   sample_reconstruction_plots='both', verbose=True)
