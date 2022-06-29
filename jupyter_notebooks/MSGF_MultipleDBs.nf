#!/usr/bin/env nextflow
// Runs MSGF+ on a set of samples against every database in the specified folder.

// Parameters supplied by the user
params.spectra_folder = false
params.databases_folder = false
params.mod = false
params.output_folder = false

modFile = file(params.mod)

spectrumFileChannel = Channel.fromPath("${params.spectra_folder}/*_dta.txt")
dbDir = file("${params.databases_folder}")

databases = []
dbDir.eachFile { item ->
	if( item.getExtension() == "fasta" ) {
		databases << file(item)
	}
}

// Run MS-GF+ on each sample-database combination
process msgfPlus {
    container "emlee/msgf:v4.0"
    cpus 16
    memory "32 GB"
    
    input:
    file spectra_file from spectrumFileChannel
	file modFile
	each file(db) from databases
	
	output:
	file "${spectra_file.name.replaceAll(/.txt/, '')}_${db.name.replaceAll(/.fasta/, '')}.mzid" into mzid_ch
	
	afterScript "rm -r *"

	script:
    """
	java -Xmx16g -jar /root/MSGFPlus.jar -s ${spectra_file} -d ${db} -tda 1 -ntt 1 -mod ${modFile} -o ${spectra_file.name.replaceAll(/.txt/, '')}_${db.name.replaceAll(/.fasta/, '')}.mzid		
	"""

}

// Convert result files from MZID to TSV and place them in publishDir
process convertMZID {
    container "emlee/msgf:v3.0"
    cpus 2
    memory "4 GB"
    publishDir "${params.output_folder}/"
    
    input:
    file mzid_file from mzid_ch
	
	output:
	file "${mzid_file.name.replaceAll(/.mzid/, '')}.tsv"
	
	afterScript "rm -r *"

	script:
    """
	java -Xmx3g -cp /root/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i ${mzid_file} -showDecoy 1
	"""

}
