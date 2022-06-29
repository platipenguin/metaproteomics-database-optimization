#!/usr/bin/env nextflow
// Runs MSGF+ on a set of samples against every database that matches the spectra's ID number.

// Parameters supplied by the user
params.spectra_folder = false
params.databases_folder = false
params.mod = false
params.output_folder = false

modFile = file(params.mod)
dbDir = file("${params.databases_folder}")
spectraDir = file("${params.spectra_folder}")

// Pair spectra and databases based on their ID numbers, then pass them to msgfPlus as a tuple
// Assumes first four characters of spectrum file name are the sample ID
pairList = []
spectraDir.eachFile { spectraFile ->
	id = spectraFile.getSimpleName()[0..3]
	dbDir.eachFile { dbFile ->
		if( dbFile.getSimpleName().indexOf(id) != -1) {
			pairList << [spectraFile, dbFile]
		}
	}
}

executionChannel = Channel.fromList(pairList)

// Run MS-GF+ on each sample-database combination
process msgfPlus {
    container "emlee/msgf:v4.0"
    cpus 16
    memory "32 GB"
    
    input:
    tuple file(spectra_file), file(db) from executionChannel
	file modFile
	
	output:
	file "${db.name.replaceAll(/.fasta/, '')}_${spectra_file.name.replaceAll(/.txt/, '')}.mzid" into mzid_ch
	
	afterScript "rm -r *"

	script:
    """
	java -Xmx16g -jar /root/MSGFPlus.jar -s ${spectra_file} -d ${db} -tda 1 -ntt 1 -mod ${modFile} -o ${db.name.replaceAll(/.fasta/, '')}_${spectra_file.name.replaceAll(/.txt/, '')}.mzid		
	"""

}

// Convert result files from MZID to TSV and place them in publishDir
process convertMZID {
    container "emlee/msgf:v4.0"
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