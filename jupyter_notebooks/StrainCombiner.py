# This program is used to compile the protein information from multiple strains of a species into a single .fasta file.
# The program ignores any files that are not .fasta files.
# Files will be named 'foldername'_combined.fasta
# ARGUMENT1 - Directory containing folders, each containing protein information for a single species
# ARGUMENT2 - Directory where resulting fasta files should be placed
import os, sys, math
from pathlib import Path

# FUNCTION DEFINITIONS
# Returns the name of the species for the given protein entry
def getSpecies(entry):
	startIndex = entry.find('OS=') + 3
	endIndex = entry.find('OX=') if entry.find('OS=') != -1 else entry.find('\n')
	return entry[startIndex:endIndex]

# MAIN PROGRAM THREAD
# Check that the user-entered arguments are valid
if len(sys.argv) < 2:
	raise Exception('You must enter a directory name')
seqDirectory = Path(sys.argv[1])
if not seqDirectory.exists():
	raise Exception(f'Could not find the specified directory {seqDirectory}')
if len(sys.argv) < 3:
	raise Exception('You must provide a directory where output files can be placed.')
outDirectory = Path(sys.argv[2])
if not outDirectory.exists():
	raise Exception(f'Could not find the specified directory {outDirectory}')

# Iterate through each folder and create a single .fasta file output
for dir in seqDirectory.iterdir():
	# Compile together all the identical sequences
	seqDict = {}
	nameDict = {}
	taxaIdentifier = dir.name
	taxaHolder = taxaIdentifier.split('_')
	taxaName = ''
	for txt in taxaHolder:
		taxaName += txt + ' '
	for f in dir.iterdir():
		if f.suffix != '.fasta':
			continue
		rawText = ''
		with f.open(mode='r') as infile:
			rawText = infile.read()
		entries = rawText.split('\n>')
		entries[0] = entries[0][1:]
		for entry in entries:
			firstSpace = entry.find(' ')
			firstNewline = entry.find('\n')
			osStart = entry.find('OS=')
			oxStart = entry.find('OX=')
			id = entry[:firstSpace]
			name = entry[firstSpace:osStart - 1]
			species = entry[osStart + 3:oxStart - 1]
			if taxaName.find('OX=') == -1:
				taxaName += entry[oxStart:firstNewline]
			seq = entry[firstNewline:].replace('\n', '')
			if not seq in seqDict.keys():
				seqDict[seq] = id
				nameDict[id] = name
	with outDirectory.joinpath(taxaIdentifier + '.fasta').open(mode='w', newline='') as newfile:
		for seq in seqDict:
			id = seqDict[seq]
			name = nameDict[id]
			newfile.write(f'>{id} {name} OS={taxaName}\n{seq}\n')
	print(f'Finished writing sequences for {taxaName}')

print('Strain Combiner completed successfully')











			
			