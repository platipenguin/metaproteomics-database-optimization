# This program downloads a list of files from an NCBI FTP site that are delineated in a csv file.
# Assumes that the CSV file has a header, and skips the first row.
# ARGUMENT1 - CSV file that contains the URLs of the files to download.
# ARGUMENT2 - The suffix of the files to download. i.e. "_genomic.fna", "_protein.faa", etc. The program will automatically append ".gz" to the url so it can get the files.
# ARGUMENT3 - (optional) Path to directory where downloaded files should be placed.
# ARGUMENT4 - (optional) 0-indexed column number in the csv file that should be used to rename each downloaded file. If none is provided, files will be named "TEMP_Genus_Species_Strain_BioSample"
# ex - python DownloadFromNCBIFTP.py penguins.csv _genomic.fna.gz Downloads\Penguins 0
import sys, os, csv, urllib.request
from pathlib import Path

# FUNCTION LIST
def getName(row, suffix):
	nameList = row[0].split(' ')
	genus = nameList[0]
	species = nameList[1]
	strain = row[2]
	sample = row[3].strip()
	return 'TEMP%' + genus + '%' + species + '%' + strain + '%' + sample + '.fasta.gz'

# Check user-entered arguments
if len(sys.argv) < 2:
	raise Exception('You must enter a csv file to parse.')
urlsFile = Path(sys.argv[1])
if not urlsFile.exists():
	raise Exception('Could not find the specified file "' + sys.argv[1])
if len(sys.argv) < 3:
	raise Exception('You must specify a file suffix.')
outPath = ''
if len(sys.argv) > 3:
	outPath = sys.argv[3] + '/'
	if not Path(sys.argv[3]).exists():
		raise Exception('Could not find output directory "' + sys.argv[3] + '"')
nameCol = -1
if len(sys.argv) > 4:
	try:
		nameCol = int(sys.argv[4])
	except:
		raise Exception('Argument 4 must be an integer')

# Get the download urls from the CSV file and determine what each download will be named
urls = {}
hasHeader = True
downloadCol = 14
with urlsFile.open() as csvFile:
	reader = csv.reader(csvFile, delimiter=',')
	for row in reader:
		if hasHeader:
			hasHeader = False
			continue
		if row[downloadCol] == '-' or row[downloadCol] == '':
			continue
		if len(row) < 2:
			break
		sampleSignifier = row[downloadCol].split('/')[len(row[downloadCol].split('/')) - 1]
		downloadPath = row[downloadCol] + '/' + sampleSignifier + sys.argv[2] + '.gz'
		fileName = row[nameCol] + '.fasta.gz' if nameCol != -1 else getName(row, sys.argv[2])
		urls[fileName] = downloadPath

# Download the files
for fileName in urls.keys():
	print('Downloading ' + urls[fileName])
	try:
		urllib.request.urlopen(urls[fileName])
		print(fileName)
	except:
		print(f'No {sys.argv[2]} file for {fileName}')
		continue
	req = urllib.request.Request(urls[fileName])
	data = None
	with urllib.request.urlopen(req) as response:
		data = response.read()
	with open(outPath + fileName, 'wb') as newFile:
		newFile.write(data)

print('Download script completed successfully.')




















