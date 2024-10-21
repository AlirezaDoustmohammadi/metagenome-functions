#!/usr/bin/env nextflow

params.biom_table = "/home/git/picrust/input/feature-table.biom"
params.seq = "/home/git/picrust/input/representative-sequences.fasta"
params.proj_dir = "/home/git/picrust/"
params.cpu = 10


log.info """\
         PICRUSt2 Pipeline    
         ===================================
         Raw Biom Table File: 	 : ${params.biom_table}
		 Raw Sequence File       : params.seq
         Output Directory        : ${params.proj_dir}/picrust2_out_pipeline
         """
         .stripIndent()
         
workflow {
		
	place_reads_into_refe_tree().set{place_reads_ch}
	
	hidden_state_prediction(place_reads_ch).set{hidden_state_prediction_ch}
	
	generate_metagenome_predictions(hidden_state_prediction_ch).set{metagenome_predictions_ch}
	
	pathway_level_inference(metagenome_predictions_ch)

}


workflow.onComplete {
    def msg = """\
        Pipeline execution summary ${workflow.scriptName}
        ---------------------------
	
	ScriptName		: ${workflow.scriptName}
	projectDir		: ${projectDir}
        workDir     : ${workflow.workDir}
		command     : ${workflow.commandLine}
        Ended at	: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
		container	: ${container}
        exit status : ${workflow.exitStatus}
		error		: ${workflow.errorMessage}
        """
        .stripIndent()
}


workflow.onError {
    def msg = """\
        Pipeline execution summary ${workflow.scriptName}
        ---------------------------
		ScriptName	: ${workflow.scriptName}
		projectDir	: ${projectDir}
        workDir     : ${workflow.workDir}
		command     : ${workflow.commandLine}
        Success     : ${workflow.success}
		container	: ${container}
        exit status : ${workflow.exitStatus}
		error		: ${workflow.errorMessage}
		errorReport	: ${workflow.errorReport}
        """
        .stripIndent()
}


process place_reads_into_refe_tree {

	conda 'picrust2'
	
	tag "Place reads into reference tree"
		
	output:
		tuple val("${params.proj_dir}/picrust2_out_pipeline/out.tre")
		
		
	script:
			"""
            mkdir -p ${params.proj_dir}/picrust2_out_pipeline
            place_seqs.py -s ${params.seq} -o ${params.proj_dir}/picrust2_out_pipeline/out.tre -p ${params.cpu}\
						  --intermediate ${params.proj_dir}/picrust2_out_pipeline/intermediate/place_seqs
			
			"""
}


process hidden_state_prediction {

	tag "Hidden-state prediction of gene families"
	
	input:
		val (ch)
		
	output:
		tuple val("${params.proj_dir}/picrust2_out_pipeline/KO_predicted.tsv.gz")
		
		
	script:
        """		
		
		
		hsp.py -i 16S -t ${params.proj_dir}/picrust2_out_pipeline/out.tre \
			   -o ${params.proj_dir}/picrust2_out_pipeline/marker_predicted_and_nsti.tsv.gz -p ${params.cpu} -n

		hsp.py -i EC -t ${params.proj_dir}/picrust2_out_pipeline/out.tre \
			   -o ${params.proj_dir}/picrust2_out_pipeline/EC_predicted.tsv.gz -p ${params.cpu}
		
		hsp.py -i KO -t ${params.proj_dir}/picrust2_out_pipeline/out.tre \
			   -o ${params.proj_dir}/picrust2_out_pipeline/KO_predicted.tsv.gz -p ${params.cpu}
		
        """
}


process generate_metagenome_predictions {

	tag "Generate metagenome predictions"

	input:
		val (ch)
		
	output:
		val ("${params.proj_dir}/picrust2_out_pipeline/KO_metagenome_out")
	
	script:
	"""
	
		  metagenome_pipeline.py -i ${params.biom_table} \
								 -m ${params.proj_dir}/picrust2_out_pipeline/marker_predicted_and_nsti.tsv.gz \
								 -f ${params.proj_dir}/picrust2_out_pipeline/EC_predicted.tsv.gz \
								 -o ${params.proj_dir}/picrust2_out_pipeline/EC_metagenome_out --strat_out
								 
								 
		  metagenome_pipeline.py -i ${params.biom_table} \
								 -m ${params.proj_dir}/picrust2_out_pipeline/marker_predicted_and_nsti.tsv.gz \
								 -f ${params.proj_dir}/picrust2_out_pipeline/KO_predicted.tsv.gz \
								 -o ${params.proj_dir}/picrust2_out_pipeline/KO_metagenome_out --strat_out
					
	"""
}


process pathway_level_inference {

	// Tag the process
	tag "Pathway-level inference"
	
	input:
		val (ch)
		
	script:
        """
		pathway_pipeline.py -i ${params.proj_dir}/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
						    -o ${params.proj_dir}/picrust2_out_pipeline/pathways_out -p ${params.cpu}
						
		# Add functional descriptions
		add_descriptions.py -i ${params.proj_dir}/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
						    -m EC \
							-o ${params.proj_dir}/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

		add_descriptions.py -i ${params.proj_dir}/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv.gz \
							-m METACYC \
							-o ${params.proj_dir}/picrust2_out_pipeline/pathways_out/path_abun_unstrat_descrip.tsv.gz
  
        """
}