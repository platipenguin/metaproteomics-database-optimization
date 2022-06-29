# This program iterates through a folder containing .fasta protein sequence files downloaded from NCBI.
# The program reads through a file, removes the '[Group]' piece in the header of each entry, and replaces it with 'OS=SpeciesName OX=SpeciesID' to make the files conform with uniprot .fasta files.
# Each file's name must follow the following naming convention: TEMP%Genus%species%Strain%IDNumber.fasta
# ARGUMENT1 - Path to folder containing files to prepare.
import os, sys
from pathlib import Path

# FUNCTION LIST
# Takes in a file name of the format 'TEMP%Genus%species%Strain%IDNumber.fasta' and returns a dictionary with keys 'genus', 'species', 'strain', and 'id'.
def extractTaxaInfo(fileName):
	vals = fileName.split('%')
	if len(vals) != 5 and len(vals) != 4:
		raise Exception('Found TEMP file with incorrect filename formatting - ' + fileName)
	id = vals[4][:vals[4].rfind('.')] if len(vals) == 5 else 'x'
	strain = vals[3] if len(vals) == 5 else vals[3][:vals[3].rfind('.')]
	return {'genus':vals[1], 'species':vals[2], 'strain':strain, 'id':id}

# MAIN PROGRAM THREAD
# Check that the user entered a valid directory name
if len(sys.argv) < 2:
	print ('Error: You must enter a results directory and a database file.')
	sys.exit()
resPath = Path(sys.argv[1])
if not resPath.exists():
	raise Exception(f'Could not find the specified folder "{sys.argv[1]}"')

# Create a new database file for each 'TEMP' species file in the folder
for file in resPath.iterdir():
	if file.name[:4] != 'TEMP': # Ignore files without 'TEMP' at the beginning of their names
		continue
	taxaInfo = extractTaxaInfo(file.name)
	replacementStr = f"OS={taxaInfo['genus']} {taxaInfo['species']} {taxaInfo['strain']} OX={taxaInfo['id']}"
	proteinList = []
	with file.open(mode='r') as infile:
		rawText = infile.read()
		proteinList = rawText.split('\n>')
		proteinList[0] = proteinList[0][1:]
		del rawText
		for i in range(len(proteinList)):
			startIndex = proteinList[i].rfind('[')
			endIndex = proteinList[i].rfind(']') + 1
			proteinList[i] = proteinList[i][:startIndex] + replacementStr + proteinList[i][endIndex:]
	with resPath.joinpath(f"{taxaInfo['genus']}{taxaInfo['species'].capitalize()}_{taxaInfo['id']}.fasta").open(mode='w', newline='') as outfile:
		for text in proteinList:
			outfile.write(f'>{text}\n')

print('Species Sequence Prepper completed successfully')