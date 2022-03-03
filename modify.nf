sdrf_file = params.input  
ch_sdrf = Channel.fromPath(sdrf_file, checkIfExists: true)                                                                          
process sdrf_parsing {

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
    file sdrf from ch_sdrf

    output:
    file "openms.tsv" into ch_sdrf_config_openms, ch_sdrf_config_file
	file "experimental_design.tsv" into ch_trans_experiment

    when:
    sdrf_file

    script:
    """
    ## -t2 since the one-table format parser is broken in OpenMS2.5
    ## -l for legacy behavior to always add sample columns
    parse_sdrf convert-openms -t2 -l -s ${sdrf} > sdrf_parsing.log
    """
}

ch_sdrf_config_file
.splitCsv(skip: 1, sep: '\t')
.multiMap{ row -> id = row.toString().md5()
                comet_settings: msgf_settings: tuple(id,
                                row[2],
                                row[3],
                                row[4],
                                row[5],
                                row[6],
                                row[7],
                                row[8],
                                row[9],
                                row[10])
                idx_settings: tuple(id,                                                                                                  
                                row[10])
                luciphor_settings:                                                                                                     
                                tuple(id,
                                row[9])
                mzmls: tuple(id, !params.root_folder ? 
                                row[0] :
                                params.root_folder + "/" + (params.local_input_type ?
                                    row[1].take(row[1].lastIndexOf('.')) + '.' + params.local_input_type :
                                    row[1]))}
.set{ch_sdrf_config}

ch_sdrf_config.mzmls
.branch {
        raw: hasExtension(it[1], 'raw')
        mzML: hasExtension(it[1], 'mzML')
}
.set {branched_input}


branched_input.mzML
.branch {
    nonIndexedMzML: file(it[1]).withReader {
                        f = it;
                        1.upto(5) {
                            if (f.readLine().contains("indexedmzML")) return false;
                        }
                        return true;
                    }
    inputIndexedMzML: file(it[1]).withReader {
                        f = it;
                        1.upto(5) {
                            if (f.readLine().contains("indexedmzML")) return true;
                        }
                        return false;
                    }
}
.set {branched_input_mzMLs}

process raw_file_conversion {

    label 'process_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, path(rawfile) from branched_input.raw

    output:
     tuple mzml_id, file("*.mzML") into mzmls_converted

    script:
     """
     ThermoRawFileParser.sh -i=${rawfile} -f=2 -o=./ > ${rawfile}_conversion.log
     """
}

process PTM_parse{
	input:
		file openms from ch_sdrf_config_openms

	output:
		file "params.csv" into Dia_ch, ch_diann_params

	script:
		"""
		python ${params.openms_parse} ${openms} ${params.unimod} > params.tsv
		"""
}
	
fasta_file = params.database  
ch_fasta = Channel.fromPath(fasta_file, checkIfExists: true)   

process Generate_spectral_library{
	input:
		file fasta from ch_fasta
	
	output:
		file 'lib.predicted.speclib' into Dia_database_search 
	
	"""
	${params.diann} \\
			--fasta ${params.database} \\
			--fasta-search \\
			--gen-spec-lib \\
			--smart-profiling \\
			--predictor
	"""
}

process Dia_NN{

	input:	
		file Params from Dia_ch
		file speclib from Dia_database_search
		file "mzMLs/*" from mzmls_converted.collect()

	output:
		file "report_*.tsv" into ch_trans_report

	script:
		
		"""
		#!/bin/bash
		for file in \$(ls ./mzMLs/)
		do
            	if [ \${file##*.} = "mzML" ];
            	then
                mkdir \${file%%.*}_mzML
                cp ./mzMLs/\${file} \${file%%.*}_mzML
                FixedPTM=`cat $Params | grep \${file%%.*} | awk -F'\t' '{print \$3}'`
                VariablePTM=`cat $Params | grep \${file%%.*} | awk -F'\t' '{print \$4}'`		
                PMTunit=`cat $Params | grep \${file%%.*} | awk -F'\t' '{print \$7}'`
                FMPunit=`cat $Params | grep \${file%%.*} | awk -F'\t' '{print \$9}'`
                Enzyme=`cat $Params | grep \${file%%.*} | awk -F'\t' '{print \$11}'`

                if [ \$Enzyme = "Trypsin" ]; 
                then
                    enzyme="K*,R*,!*P"
                fi
                if [ \$Enzyme = "Arg-C" ]; 
                then
                    enzyme="R*,!*P"
                fi
                if [ \$Enzyme = "Asp-N" ]; 
                then
                    enzyme="*B,*D"
                fi
                if [ \$Enzyme = "Chymotrypsin" ]; 
                then
                    enzyme="F*,W*,Y*,L*,!*P"
                fi
                if [ \$Enzyme = "Lys-C" ]; 
                then
                    enzyme="K*,!*P"
                fi

                if [ \$PMTunit = "ppm" ]; 
                then
                    PMT=`cat $Params | grep \${file%%.*} | awk -F'\t' '{print \$6}'`
                else
                    PMT=${params.precursor_mass_tolerance}
                fi

                if [ \$FMPunit = "ppm" ]; 
                then
                    FMP=`cat $Params | grep \${file%%.*} | awk -F'\t' '{print \$8}'`
                else
                    FMP=${params.fragment_mass_tolerance}
                fi

                ${params.diann} \\
                    --dir \${file%%.*}_mzML \\
                    --lib ./lib.predicted.speclib \\
                    --out report_\${file%%.*}.tsv \\
                    --no-stats \\
                    --mass-acc-ms1 \${PMT} \\
                    --mass-acc \${FMP} \\
                    --cut \${enzyme} \\
                    --clear-mods \\
                    --full-unimod \\
                    --smart-profiling \\
                    --reanalyse \\
                    \${FixedPTM} \\
                    \${VariablePTM} \\
                    --missed-cleavages ${params.allowed_missed_cleavages} \\
                    --min-pep-len ${params.min_peptide_length} \\
                    --max-pep-len ${params.max_peptide_length} \\
                    --min-pr-charge ${params.min_precursor_charge} \\
                    --max-pr-charge ${params.max_precursor_charge} \\
                    --min-pr-mz ${params.min_pr_mz} \\
                    --max-pr-mz ${params.max_pr_mz} \\
                    --min-fr-mz ${params.min_pr_mz} \\
                    --max-fr-mz ${params.max_pr_mz} \\
                    --threads ${params.max_cpus} \\
                    --var-mods ${params.max_var_mods}
            	fi
		done
		"""						
}

process Trans_to_MSstats{
	input:
	file expriment from ch_trans_experiment
	file "Reps/*.tsv" from ch_trans_report
	"""
	#!/bin/bash  
	
	for file in \$(ls ./Reps/);  
	do
	sed '1d' ./Reps/\${file} >> ./report.tsv;
	done

	content=`cat ./report.tsv`;
	echo -e "File.Name\tRun\tProtein.Group\tProtein.Ids\tProtein.Names\tGenes\tPG.Quantity\tPG.Normalised\tPG.MaxLFQ\tGenes.Quantity\tGenes.Normalised\tGenes.MaxLFQ\tGenes.MaxLFQ.Unique\tModified.Sequence\tStripped.Sequence\tvPrecursor.Id\tPrecursor.Charge\tQ.Value\tPEP\tGlobal.Q.Value\tProtein.Q.Value\tPG.Q.Value\tGlobal.PG.Q.Value\tGG.Q.Value\tTranslated.Q.Value\tProteotypic\tPrecursor.Quantity\tPrecursor.Normalised\tPrecursor.Translated\tMs1.Translated\tQuantity.Quality\tRT\tRT.Start\tRT.Stop\tiRT\tPredicted.RT\tPredicted.iRT\tFirst.Protein.Description\tLib.Q.Value\tLib.PG.Q.Value\tMs1.Profile.Corr\tMs1.Area\tEvidence\tSpectrum.Similarity\tMass.Evidence\tCScore\tDecoy.Evidence\tDecoy.CScore\tFragment.Quant.Raw\tFragment.Quant.Corrected\tFragment.Correlations\tMS2.Scan\tIM\tiIM\tPredicted.IM\tPredicted.iIM\n\$content" > ./report.tsv

	python $params.trans_to_mastats report.tsv $params.unimod experimental_design.tsv
	"""
}
//--------------------------------------------------------------- //
//---------------------- Utility functions  --------------------- //
//--------------------------------------------------------------- //

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Check class of an Object for "List" type
boolean isCollectionOrArray(object) {
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}
