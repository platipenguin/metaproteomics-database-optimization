#!/bin/bash

set -e

module load nextflow

export NXF_VER=20.07.1

WORKING_DIR=""

nextflow \
	-C ~/nextflow.config \
	run \
	MSGF_TailoredDBs.nf \
		--spectra_folder "SampleData/R2/" \
		--databases_folder "12-16-21_NextflowMSGF_Tailored4_Combined/databases/" \
		--mod "MSGFDB_PartTryp_MetOx_20ppmParTol.txt" \
		--output_folder "1-19-21_NextflowMSGF_R2_Tailored/output/" \
		-w "$WORKING_DIR" \
		-process.queue "default" \
		-with-report 1-18-21_NextflowMSGF_R1_Tailored/Tailored_R1_Report.html \
		-resume
