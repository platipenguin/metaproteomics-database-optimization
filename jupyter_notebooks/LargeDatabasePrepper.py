# This program is used to create a fasta database file from a folder of protein sequences.
# It is specially designed for databases >250MB, which is the rough upper-bound for an MS-GF+ database.
# It reads through and copies every sequence in each .fasta file in the folder, copying them out to the new database file.
# However, each time the output file reaches roughly exceeds 250MB, the program starts writing to a new file, making a sequence of database files, i.e. db-1.fasta, db-2.fasta, etc.
# The program DOES NOT verify that there are no proteins with identical sequences in the databases. Data going in must be non-redundant.
# ARGUMENT1 - Path to folder with sequences to place in the database.
# ARGUMENT2 - What the new database should be named (do not include the ".fasta" suffix).
# ARGUMENT3 - Path to directory where the database files should be written (must be different than the path to the sequences folder).
# ex - python LargeDatabasePrepper.py Sequences\DBBuilder MyPreppedDatabase Sequences\CompiledDatabases
import sys
from pathlib import Path

# FUNCTION LIST
# Object used to store information about a protein
class Protein:
	def __init__(self, entry):
		self.id = self.extractProteinID(entry)
		self.isContaminant = self.id.find('Contaminant_') != -1
		self.taxa = self.extractSpecies(entry)
		self.name = self.extractProteinName(entry)
		self.data = 'OX=x' if self.isContaminant else entry[entry.find('OX='):entry.find('\n')]
		self.sequence = entry[entry.find('\n') + 1:]
		self.hitPeptides = set()
	
	# This function takes a single protein entry from a .fasta file and returns all the species for that entry
	def extractSpecies(self, entry):
		if self.isContaminant:
			return 'Contaminant'
		stringStart = entry.find('OS=') + 3
		stringEnd = entry.find(' OX=')
		return entry[stringStart:stringEnd].split(';')

	# This function takes a single protein entry from a .fasta file and returns the protein identifier for that entry.
	def extractProteinID(self, entry):
		stringEnd = entry.find(' ')
		return entry[0:stringEnd]

	# This function takes a single protein entry from a .fasta file and returns the protein identifier for that entry.
	def extractProteinName(self, entry):
		stringStart = entry.find(' ') + 1
		stringEnd = entry.find('(') if self.isContaminant else entry.find(' OS=')
		return entry[stringStart:stringEnd]
	
	# Get the fasta entry for this protein, including all taxa that encode it
	def getEntry(self):
		taxaEntry = self.taxa[0]
		if len(self.taxa) > 1:
			for i in range(1, len(self.taxa)):
				taxaEntry = f'{taxaEntry};{self.taxa[i]}'
		return f'>{self.id} {self.name} OS={taxaEntry} {self.data}\n{self.sequence}\n'	

# MAIN PROGRAM THREAD
# Check that the user entered a valid directory name
if len(sys.argv) < 2:
	raise Exception('You must enter a folder name containing fasta files to compile into a database.')
sequencePath = Path(sys.argv[1])
if sequencePath.exists() == False:
	raise Exception('Could not find the specified folder "' ,sys.argv[1],'"')
if len(sys.argv) < 3:
	raise Exception('You must specify a name for the new database.')
dbBaseName = sys.argv[2]
if len(sys.argv) < 4:
	raise Exception('You must specify a folder to place the new databases.')
outPath = Path(sys.argv[3])
if outPath.exists() == False:
	raise Exception('Could not find the specified folder "' ,sys.argv[3],'"')
if outPath == sequencePath:
	raise Exception('Sequence folder and output folder cannot be the same.')


# Gather sequences from all files and write them out to the new database
dbCount = 1
totalSeqs = 0
for f in [x for x in sequencePath.iterdir() if x.suffix == '.fasta']:
	toWrite = []
	with f.open(mode='r') as database:
		entry = []
		for line in database:
			if line[0] == '>' and entry != []:
				entry[0] = entry[0].replace('>', '')
				needEdit = entry[0].find('[') != -1
				if needEdit:
					entry[0] = entry[0].replace('[[', 'OS=')
					entry[0] = entry[0].replace('[', 'OS=')
					entry[0] = entry[0].replace(']', '')
					entry[0] = entry[0].replace('\n', ' OX=x\n')
				newProt = Protein(''.join(entry))
				toWrite.append(newProt.getEntry())
				entry = [line]
			else:
				entry.append(line)
		entry[0] = entry[0].replace('>', '')
		needEdit = entry[0].find('[') != -1
		if needEdit:
			entry[0] = entry[0].replace('[[', 'OS=')
			entry[0] = entry[0].replace('[', 'OS=')
			entry[0] = entry[0].replace(']', '')
			entry[0] = entry[0].replace('\n', ' OX=x')
		newProt = Protein(''.join(entry))
		toWrite.append(newProt.getEntry())
	totalSeqs += len(toWrite)
	print(f'Done reading {f.name}')
	# Write the proteins out to file
	_iter = 0
	while _iter < len(toWrite):
		dbPath = outPath.joinpath(f'{dbBaseName}_{str(dbCount)}.fasta')
		with dbPath.open(mode='a', newline='') as outfile:
			for i in range(_iter, len(toWrite)):
				outfile.write(toWrite[i].rstrip() + '\n')
				_iter += 1
				if dbPath.stat().st_size > 300000000:
					dbCount += 1
					break
	print(f'Done reading {f.name}. {str(len(toWrite))} sequences written.')

print(f'Database Prepper completed successfully. Database contains {str(totalSeqs)} protein sequences.')