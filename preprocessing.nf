#!/usr/bin/env nextflow

params.biom_table = "/home/git/picrust/input/feature-table.biom"
params.seq = "/home/git/picrust/input/representative-sequences.fasta"
params.input_dir = "/home/git/picrust/input/"
params.proj_dir = "/home/git/picrust/"
params.min_feature_abundance = 15
params.min_samples_abundance= 1500


log.info """\
         Qiime2 - Preprocessing    
         ===================================
         Raw Biom Table File: 	 : ${params.biom_table}
		 Raw Sequence File       : params.seq
         Output Directory        : ${params.proj_dir}/preprocessed
         """
         .stripIndent()
         
workflow {
    
	biom_table_baseName = file(params.biom_table).baseName
    seq_baseName = file(params.seq).baseName
	
	min_feature_abundance = params.min_feature_abundance.toInteger()
	min_samples_abundance = params.min_samples_abundance.toInteger()
		

	convertFormat(biom_table_baseName, seq_baseName).set{res_ch}
	
	filterFeaturesAndSamplesBasedonTotalFrequency(res_ch, biom_table_baseName, 
												  min_feature_abundance, min_samples_abundance).set{res2_ch}
	
	filterSequencesAndExporttoFormats(res2_ch, seq_baseName, biom_table_baseName).set{res3_ch}
	
	convertFilteredBiomtoTsvandSummarize(res3_ch, biom_table_baseName)

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


process convertFormat {

	conda 'qiime2-metagenome-2024.5'
	
	tag "Convert Biom to tsv & Convert fastq to qza format"
	
	input:
		val (biom_table_baseName)
		val (seq_baseName)
		
	output:
		tuple val("{params.input_dir}/${biom_table_baseName}.qza"), val("{params.input_dir}/${seq_baseName}.qza")
		
		
	script:
			"""
            
            biom convert -i ${params.biom_table} -o ${params.input_dir}/${biom_table_baseName}.tsv \
                 --to-tsv --header-key taxonomy
            # summariese the biome file
            biom summarize-table -i ${params.biom_table} -o ${params.input_dir}/${biom_table_baseName}.summary
            python ${params.proj_dir}/biomSummary.py --biom_file ${params.biom_table}
            
            
            
            # Convert feature tables from biom files to QZA-format for filtering feature table using qiime
            qiime tools import --input-path ${params.biom_table} --type 'FeatureTable[Frequency]' \
                        --input-format BIOMV210Format --output-path ${params.input_dir}/${biom_table_baseName}.qza
                        
            # Convert ASVs sequence file from fasta files to QZA-format for filtering them using qiime
            qiime tools import --type 'FeatureData[Sequence]' --input-path ${params.seq} \
                               --output-path ${params.input_dir}/${seq_baseName}.qza
			
			"""
}


process filterFeaturesAndSamplesBasedonTotalFrequency {

	tag "Filter Features And Samples Based on Total Frequency"
	
	input:
		val (res_ch)
		val (biom_table_baseName)
		val (min_feature_abundance)
		val (min_samples_abundance)
		
	output:
		tuple val("${params.proj_dir}/preprocessing/${biom_table_baseName}.filtered_features_based_total.frequency.qza")
		
		
	script:
        """		
		mkdir -p ${params.proj_dir}/preprocessing/
		
		# remove all features with a total abundance (summed across all samples) of less than threshold
		qiime feature-table filter-features --i-table ${params.input_dir}/${biom_table_baseName}.qza \
				--p-min-frequency ${min_feature_abundance} \
				--o-filtered-table \
				${params.proj_dir}/preprocessing/${biom_table_baseName}.filtered_features_based_total.frequency.qza
	
		
		# filtering samples whose total frequency is less than threshold
		qiime feature-table filter-samples \
				--i-table \
					${params.proj_dir}/preprocessing/${biom_table_baseName}.filtered_features_based_total.frequency.qza \
				--p-min-frequency ${min_samples_abundance} \
				--o-filtered-table \
					${params.proj_dir}/preprocessing/${biom_table_baseName}.filtered_samples_based_total.frequency.qza
        """
}


process filterSequencesAndExporttoFormats {

	tag "Filter Sequences And Export to Formats"

	input:
		val (res2_ch)
		val (seq_baseName)
		val (biom_table_baseName)
		
	output:
		val ("${params.proj_dir}/preprocessing/${biom_table_baseName}.filtered_samples_based_total.frequency.biom")
	
	script:
	"""
	
		  # filter sequences to match the filtered feature table 
		  qiime feature-table filter-seqs \
			--i-data ${params.input_dir}/${seq_baseName}.qza \
			--i-table ${params.proj_dir}/preprocessing/${biom_table_baseName}.filtered_samples_based_total.frequency.qza \
			--o-filtered-data ${params.proj_dir}/preprocessing/${seq_baseName}.filtered.qza

		  # convert filtered feature table & ASVs sequence to biom & fasta format
		  # convert filtered ASVs sequence file to fasta format
		  qiime tools export --input-path ${params.proj_dir}/preprocessing/${seq_baseName}.filtered.qza \
			--output-path ${params.proj_dir}/preprocessing/data_pre_processing
		
	
		  mv ${params.proj_dir}/preprocessing/data_pre_processing/dna-sequences.fasta \
				${params.proj_dir}/preprocessing/${seq_baseName}.filtered.fasta

		  # convert feature table to biom format
		  qiime tools export \
			--input-path ${params.proj_dir}/preprocessing/${biom_table_baseName}.filtered_samples_based_total.frequency.qza \
			--output-path ${params.proj_dir}/preprocessing/${biom_table_baseName}.filtered_samples_based_total.frequency.biom \
			--output-format BIOMV210Format	
					
	"""
}


process convertFilteredBiomtoTsvandSummarize {

	// Tag the process
	tag "convert Filtered Biom to Tsv and Summarize"

	// Define the input
	input:
		val (res3_ch)
		val (biom_table_baseName)
	
	script:
        """
		mkdir -p ${params.proj_dir}/preprocessing/summary/
		
		# convert filtered biom file to tsv
		biom convert -i ${params.proj_dir}/preprocessing/${biom_table_baseName}.filtered_samples_based_total.frequency.biom\
					-o ${params.proj_dir}/preprocessing/summary/${biom_table_baseName}.filtered_samples_based_total.frequency.tsv \
						--to-tsv --header-key taxonomy
						
		# summariese the biome file
		biom summarize-table -i \
			${params.proj_dir}/preprocessing/${biom_table_baseName}.filtered_samples_based_total.frequency.biom \
			-o ${params.proj_dir}/preprocessing/summary/${biom_table_baseName}.filtered_samples_based_total.frequency.summary
			
		python ${params.proj_dir}/biomSummary.py --biom_file \
				${params.proj_dir}/preprocessing/${biom_table_baseName}.filtered_samples_based_total.frequency.biom
  
        """
}