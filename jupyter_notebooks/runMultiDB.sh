#!/bin/bash

set -e

module load nextflow

WORKING_DIR=""

nextflow \
	-C ~/nextflow.config \
	run \
	MSGF_MultipleDBs.nf \
		--spectra_folder "SampleData/Subset/" \
		--databases_folder "12-21-21_NextflowMSGF_Community5_Combined/databases_subset_refined/" \
		--mod "MSGFDB_PartTryp_MetOx_20ppmParTol.txt" \
		--output_folder "12-21-21_NextflowMSGF_Community5_Combined/output_subset_refined/" \
		-w "$WORKING_DIR" \
		-process.queue "default" \
		-with-report 12-21-21_NextflowMSGF_Community5_Combined/Community5_Subset_Refined_Report.html \
		-resume
