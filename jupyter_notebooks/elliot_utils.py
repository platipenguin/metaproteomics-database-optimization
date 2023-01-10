# This file defines a number of functions and classes that are useful for parsing through and analyzing proteomics data.
import os, sys, csv, math, shutil, subprocess, multiprocessing.pool, functools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from pathlib import Path
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
sns.set(font_scale=1.3)
sns.set_style('ticks')
sns.set_palette('bright')
csv.field_size_limit(2147483647)

# Tuple containing strings with the number of each sample, in order
SAMPLE_NAMES = ('1284', '1289', '1290', '1292', '1294', '1296', '1299', '1303', '1304', '1306', '1310', '1314', '1316', '1318', '1320', '1322', '1324', '1326', '1328', '1334', '1336', '1338', '1340', '1342', '1346', '1348', '1350', '1356', '1358')
# Tuple corresponding to the BV status of each sample, in order
BV_STATUS = ('-', '-', '-', '+', '+', '+', '+', '+', '+', '+', '+', '-', '-', '+', '+', '+', '+', '-', '+', '-', '+', '-', '-', '+', '+', '+', '+', '+', '+')
SAMPLE_MANUSCRIPT_IDS = ('BV-_1','BV-_2','BV-_3','BV+_1','BV+_2','BV+_3','BV+_4','BV+_5','BV+_6','BV+_7','BV+_8','BV-_4','BV-_5','BV+_9','BV+_10','BV+_11','BV+_12','BV-_6','BV+_13','BV-_7','BV+_14','BV-_8','BV-_9','BV+_15','BV+_16','BV+_17','BV+_18','BV+_19','BV+_20')
# Tuple containing total qPCR counts for each sample in order
BACTERIAL_LOADS = (897000000,25200000,1015500000,28650000000,22050000000,8010000000,6345000000,3015000000,1500000000,4320000000,604500000,3165000000,810000000,22950000000,5190000000,2130000000,2835000000,1318500000,189000000,750000000,1060500000,39000000,217500000,1219500000,1137000000,4155000000,1041000000,540000000,751500000)
# Tuple containing a very long list of hex codes for distinct colors
COLORS = ('#FFFF00', '#1CE6FF', '#FF34FF', '#FF4A46', '#008941', '#006FA6', '#A30059', '#FFDBE5', '#7A4900', '#0000A6', '#63FFAC', '#B79762', '#004D43', '#8FB0FF', '#997D87','#5A0007', '#809693', '#FEFFE6', '#1B4400', '#4FC601', '#3B5DFF', '#4A3B53', '#FF2F80', '#61615A', '#BA0900', '#6B7900', '#00C2A0', '#FFAA92', '#FF90C9', '#B903AA', '#D16100','#DDEFFF', '#000035', '#7B4F4B', '#A1C299', '#300018', '#0AA6D8', '#013349', '#00846F', '#372101', '#FFB500', '#C2FFED', '#A079BF', '#CC0744', '#C0B9B2', '#C2FF99', '#001E09', '#00489C', '#6F0062', '#0CBD66', '#EEC3FF', '#456D75', '#B77B68', '#7A87A1', '#788D66','#885578', '#FAD09F', '#FF8A9A', '#D157A0', '#BEC459', '#456648', '#0086ED', '#886F4C','#34362D', '#B4A8BD', '#00A6AA', '#452C2C', '#636375', '#A3C8C9', '#FF913F', '#938A81','#575329', '#00FECF', '#B05B6F', '#8CD0FF', '#3B9700', '#04F757', '#C8A1A1', '#1E6E00','#7900D7', '#A77500', '#6367A9', '#A05837', '#6B002C', '#772600', '#D790FF', '#9B9700','#549E79', '#FFF69F', '#201625', '#72418F', '#BC23FF', '#99ADC0', '#3A2465', '#922329','#5B4534', '#FDE8DC', '#404E55', '#0089A3', '#CB7E98', '#A4E804', '#324E72', '#6A3A4C','#83AB58', '#001C1E', '#D1F7CE', '#004B28', '#C8D0F6', '#A3A489', '#806C66', '#222800','#BF5650', '#E83000', '#66796D', '#DA007C', '#FF1A59', '#8ADBB4', '#1E0200', '#5B4E51','#C895C5', '#320033', '#FF6832', '#66E1D3', '#CFCDAC', '#D0AC94', '#7ED379', '#012C58','#7A7BFF', '#D68E01', '#353339', '#78AFA1', '#FEB2C6', '#75797C', '#837393', '#943A4D','#B5F4FF', '#D2DCD5', '#9556BD', '#6A714A', '#001325', '#02525F', '#0AA3F7', '#E98176','#DBD5DD', '#5EBCD1', '#3D4F44', '#7E6405', '#02684E', '#962B75', '#8D8546', '#9695C5','#E773CE', '#D86A78', '#3E89BE', '#CA834E', '#518A87', '#5B113C', '#55813B', '#E704C4','#00005F', '#A97399', '#4B8160', '#59738A', '#FF5DA7', '#F7C9BF', '#643127', '#513A01','#6B94AA', '#51A058', '#A45B02', '#1D1702', '#E20027', '#E7AB63', '#4C6001', '#9C6966','#64547B', '#97979E', '#006A66', '#391406', '#F4D749', '#0045D2', '#006C31', '#DDB6D0','#7C6571', '#9FB2A4', '#00D891', '#15A08A', '#BC65E9', '#FFFFFE', '#C6DC99', '#203B3C','#671190', '#6B3A64', '#F5E1FF', '#FFA0F2', '#CCAA35', '#374527', '#8BB400', '#797868','#C6005A', '#3B000A', '#C86240', '#29607C', '#402334', '#7D5A44', '#CCB87C', '#B88183','#AA5199', '#B5D6C3', '#A38469', '#9F94F0', '#A74571', '#B894A6', '#71BB8C', '#00B433','#789EC9', '#6D80BA', '#953F00', '#5EFF03', '#E4FFFC', '#1BE177', '#BCB1E5', '#76912F','#003109', '#0060CD', '#D20096', '#895563', '#29201D', '#5B3213', '#A76F42', '#89412E','#1A3A2A', '#494B5A', '#A88C85', '#F4ABAA', '#A3F3AB', '#00C6C8', '#EA8B66', '#958A9F','#BDC9D2', '#9FA064', '#BE4700', '#658188', '#83A485', '#453C23', '#47675D', '#3A3F00','#061203', '#DFFB71', '#868E7E', '#98D058', '#6C8F7D', '#D7BFC2', '#3C3E6E', '#D83D66','#2F5D9B', '#6C5E46', '#D25B88', '#5B656C', '#00B57F', '#545C46', '#866097', '#365D25','#252F99', '#00CCFF', '#674E60', '#FC009C', '#92896B')

# Constant values for row indices where useful information can be found in .tsv result files
SPEC_ID = 1
PEPTIDE = 8
PROTEIN_HITS = 9
MSGF_SCORE = 11
SPEC_PROBABILITY = 12
Q_VALUE = 14

# Constant path to figures folder
FIGURES = Path('figures/')

# Constant paths to final result and refined database directories
GLOBAL_SUBSET_RESULTS = Path.cwd().joinpath('../8-4-20_NextflowMSGF_Combined_Global_CommunityV4_Refined/global_refined_output_condensed/')
GLOBAL_SUBSET_DB = Path.cwd().joinpath('../8-4-20_NextflowMSGF_Combined_Global_CommunityV4_Refined/global_refined_databases/')
COMMUNITY_SUBSET_RESULTS = Path.cwd().joinpath('../12-21-21_NextflowMSGF_Community5_Combined/output_subset_refined_processed/')
COMMUNITY_SUBSET_DB = Path.cwd().joinpath('../12-21-21_NextflowMSGF_Community5_Combined/databases_subset_refined/')
COMMUNITY_RESULTS = Path.cwd().joinpath('../12-21-21_NextflowMSGF_Community5_Combined/output_refined_processed/')
COMMUNITY_DB = Path.cwd().joinpath('../12-21-21_NextflowMSGF_Community5_Combined/refined_db_singlefile/')
TAILORED_SUBSET_RESULTS = Path.cwd().joinpath('../12-16-21_NextflowMSGF_Tailored4_Combined/output_subset/')
TAILORED_SUBSET_DB = Path.cwd().joinpath('../12-16-21_NextflowMSGF_Tailored4_Combined/databases_subset_refined/')
TAILORED_RESULTS = Path.cwd().joinpath('../12-16-21_NextflowMSGF_Tailored4_Combined/output_refined_processed/')
TAILORED_DB = Path.cwd().joinpath('../12-16-21_NextflowMSGF_Tailored4_Combined/databases_refined/')
POOLED_SUBSET_RESULTS = Path.cwd().joinpath('../4-21-21_NextflowMSGF_Subset_Pooled/output_refined_processed/')
POOLED_SUBSET_DB = Path.cwd().joinpath('../4-21-21_NextflowMSGF_Subset_Pooled/databases_refined/')
POOLED_RESULTS = Path.cwd().joinpath('../4-9-21_NextflowMSGF_Combined_Pooled/output_combined_refined_processed/')
POOLED_DB = Path.cwd().joinpath('../4-9-21_NextflowMSGF_Combined_Pooled/database_refined/')
SINGLE_SUBSET_RESULTS = Path.cwd().joinpath('../4-2-21_NextflowMSGF_Combined_Individual/output_subset/')
SINGLE_SUBSET_DB = Path.cwd().joinpath('../4-2-21_NextflowMSGF_Combined_Individual/databases_subset_refined/')
SINGLE_RESULTS = Path.cwd().joinpath('../4-2-21_NextflowMSGF_Combined_Individual/output_refined_processed/')
SINGLE_DB = Path.cwd().joinpath('../4-2-21_NextflowMSGF_Combined_Individual/databases_refined/')
HYBRID_SUBSET_RESULTS = Path.cwd().joinpath('../12-17-21_NextflowMSGF_Combined_Hybrid2/output_subset/')
HYBRID_SUBSET_DB = Path.cwd().joinpath('../12-17-21_NextflowMSGF_Combined_Hybrid2/databases_subset_refined/')
HYBRID_RESULTS = Path.cwd().joinpath('../12-17-21_NextflowMSGF_Combined_Hybrid2/output_refined_processed/')
HYBRID_DB = Path.cwd().joinpath('../12-17-21_NextflowMSGF_Combined_Hybrid2/databases_refined/')

# Sets of protein IDs to check the type of a specific protein
HUMAN_PROTEINS = set()
CONTAMINANT_PROTEINS = set()
BACTERIAL_PROTEINS = set()
FUNGAL_PROTEINS = set()
TRICHOMONAS_PROTEINS = set()

# Object used to store information about a protein
class Protein:
    def __init__(self, entry):
        if entry == '':
            self.initEmpty()
            return
        self.id = self.extractProteinID(entry)
        self.isContaminant = self.id.find('Contaminant_') != -1
        self.taxa = self.extractSpecies(entry)
        self.name = self.extractProteinName(entry)
        self.data = '' if self.isContaminant else entry[entry.find('OS='):entry.find('\n')]
        self.sequence = entry[entry.find('\n') + 1:].replace('\n', '')
        self.hitPeptides = set()
    
    def initEmpty(self):
        self.id = None
        self.isContaminant = False
        self.taxa = None
        self.name = None
        self.data = 'OX=x'
        self.sequence = None
        self.hitPeptides = set()
    
    # This function takes a single protein entry from a .fasta file and returns all the species for that entry
    def extractSpecies(self, entry):
        if self.isContaminant:
            return ['Contaminant']
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
        stringEnd = entry.find('\n') if self.isContaminant else entry.find(' OS=')
        return entry[stringStart:stringEnd]
    
    # Adds the specified taxa name to this protein's taxa list
    def addTaxa(self, taxa):
        if not taxa in self.taxa:
            self.taxa.append(taxa)
    
    # Get the fasta entry for this protein, including all taxa that encode it
    def getEntry(self):
        if self.isContaminant:
            return f'>{self.id} {self.name}\n{self.sequence}\n'
        taxaEntry = self.taxa[0]
        if len(self.taxa) > 1:
            for i in range(1, len(self.taxa)):
                taxaEntry = f'{taxaEntry};{self.taxa[i]}'
        return f'>{self.id} {self.name} OS={taxaEntry} {self.data}\n{self.sequence}\n'
    
    # Get the fasta entry for this protein with amino acid sequence broken into 60 character lines
    def getFormattedEntry(self):
        if self.isContaminant:
            raise Exception('Cannot get formatted entry for contaminant protein.')
        taxaEntry = self.taxa[0]
        if len(self.taxa) > 1:
            for i in range(1, len(self.taxa)):
                taxaEntry = f'{taxaEntry};{self.taxa[i]}'
        sequenceEntry = []
        for i in range(0, len(self.sequence), 60):
            if i + 60 >= len(self.sequence):
                sequenceEntry.append(self.sequence[i:] + '\n')
            else:
                sequenceEntry.append(self.sequence[i:i + 60] + '\n')
        return f'>{self.id} {self.name} OS={taxaEntry} {self.data}\n{"".join(sequenceEntry)}'
    
    def equals(self, otherProtein):
        return self.id == otherProtein.id
    
    # Returns a new Protein object that is identical to this one
    def copy(self):
        return Protein(self.getEntry()[1:])

# Processes through all the proteins in a fastafile and adds their IDs to the specified set
def fastaToIDSet(fastafile, outputset):
    data = ''
    with fastafile.open(mode='r') as infile:
        data = infile.read()
    dataList = data.split('\n>')
    dataList[0] = dataList[0][1:]
    del data
    for sequence in dataList:
        outputset.add(Protein(sequence).id)

fastaToIDSet(Path.cwd().joinpath('../PublicSequences/SwissProt_Human_12-21.fasta'), HUMAN_PROTEINS)
fastaToIDSet(Path.cwd().joinpath('../PublicSequences/Contaminant_12-21.fasta'), CONTAMINANT_PROTEINS)
fastaToIDSet(Path.cwd().joinpath('../PublicSequences/NCBI_Bacteria_12-21.fasta'), BACTERIAL_PROTEINS)
fastaToIDSet(Path.cwd().joinpath('../ShotgunMetagenomics/genes_annotated.fasta'), BACTERIAL_PROTEINS)
fastaToIDSet(Path.cwd().joinpath('../PublicSequences/NCBI_VaginalFungi_12-21.fasta'), FUNGAL_PROTEINS)
fastaToIDSet(Path.cwd().joinpath('../PublicSequences/Combined_AllNCBI_Eukaryotes_12-21/Trichomonas_vaginalis.fasta'), TRICHOMONAS_PROTEINS)

# Adds the bacterial and fungal proteins in the global database to the BACTERIAL_PROTEINS and FUNGAL_PROTEINS sets
# Should only be called when doing analysis with the Global database, unnecessary otherwise.
def initWithGlobal():
    fastaToIDSet(Path.cwd().joinpath('../PublicSequences/global_bacteria_proteins.fasta'), BACTERIAL_PROTEINS)
    fastaToIDSet(Path.cwd().joinpath('../PublicSequences/global_fungi_proteins.fasta'), FUNGAL_PROTEINS)

# Adds the bacterial proteins from all the individual bacterial fasta files to the BACTERIAL_PROTEINS sets
# Should only be called when doing analysis with the 16S_Reference database, unnecessary otherwise.
def initWithReference():
    for p in Path('C:/Users/emlee/Documents/MSGFp/Sequences/Bacteria/AllNCBI_12-19/').rglob('*'):
        if not p.suffix == '.fasta':
            continue
        with open(p, 'r') as database:
            rawText = database.read()
            dataList = rawText.split('\n>')
            dataList[0] = dataList[0][1:]
            del rawText
            for sequence in dataList:
                newProt = Protein(sequence)
                BACTERIAL_PROTEINS.add(newProt.id)

# Object representing the genome of a species that tracks what sample peptides it covers.
class Genome:
    def __init__(self, name):
        self.name = name
        self.peptideSet = set()
    
    # If the peptide isn't already in this genome's list of covered peptides, adds it to the list.
    def addPeptide(self, newPep):
        self.peptideSet.add(newPep)
    
    @property
    def numPeptides(self):
        return len(self.peptideSet)

# Object used to store information about a single spectrum hit in a result file
class Spectrum:
    def __init__(self, resRow):
        self.id = resRow[SPEC_ID]
        self.peptide = resRow[PEPTIDE]
        self.hits = getProteinHitList(resRow, 'all')
        self.q = float(resRow[Q_VALUE])
    
    def equals(self, otherSpectrum):
        return self.id == otherSpectrum.id

# Object that makes it easier to work with the data in a TSV result file. Defines rows by their spectrum ID.
class TSVFile:
    def __init__(self, inputFile=None):
        self.rows = {} # key=specID, value=list representing tsv row
        self.header = []
        if inputFile == None:
            return
        with inputFile.open(mode='r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            isFirst = True
            for row in reader:
                if isFirst:
                    self.header = row
                    isFirst = False
                else:
                    self.rows[row[SPEC_ID]] = row
    
    def length(self):
        return len(self.rows)
    
    # Get the row represented by the given specID
    def getRow(self, specID):
        return self.rows[specID]
    
    # Sets the row value for this TSV file at the specified specID
    def setRow(self, specID, newRow):
        if newRow[SPEC_ID] != specID:
            raise Exception(f'ID of new row {newRow[SPEC_ID]} does not match the specified ID {specID}')
        self.rows[specID] = newRow
    
    # Get the MSGF score as a float for the row with the specified specID
    def getMSGFScore(self, specID):
        return float(self.rows[specID][MSGF_SCORE])
    
    # Writes this object out to a TSV file.
    def writeToFile(self, outputPath):
        toWrite = []
        for r in self.rows.values():
            toWrite.append(r)
        toWrite.sort(key=lambda x: float(x[Q_VALUE]), reverse=False)
        with outputPath.open(mode='w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            writer.writerow(self.header)
            for rowToWrite in toWrite:
                writer.writerow(rowToWrite)
            

# Object used to build and access a protein database
# Excludes contaminant proteins, but includes all other types of proteins in the database.
# Uses a hash set as a data backend so that finding a protein by ID is O(1).
class ProtRef:
    # Create a new ProtRef object using the specified dbFile. If disallowRepeats = True, the constructor will throw an exception if there are proteins with identical IDs.
    def __init__(self, dbFile, disallowRepeats=True):
        self.db = {} # key=protein id, value=protein object
        self.peps = {} # key=peptide sequence, value=set of protein id's
        self.taxa = set()
        with open(dbFile, 'r') as database:
            rawText = database.read()
            dataList = rawText.split('\n>')
            dataList[0] = dataList[0][1:]
            del rawText
            for sequence in dataList:
                if sequence.find('Contaminant_') != -1:
                    continue
                newProt = Protein(sequence)
                if disallowRepeats and newProt.id in self.db.keys():
                    raise Exception('There are at least two proteins ' + str(dbFile) + ' with ID "' + newProt.id + '". Protein IDs must be unique.')
                self.addProt(newProt)
                for taxa in newProt.taxa:
                    self.taxa.add(taxa)
    
    # Adds a protein to the reference database
    def addProt(self, newProt):
        self.db[newProt.id] = newProt
    
    # Documents which proteins in this ProtRef hit which peptides in the results
    # Can take a path representing a single file, or a path representing a directory full of files.
    def setPeptideHits(self, results):
        self.peps = {}
        resType = type(results)
        if resType == list:
            for f in results:
                self._setPeptideHits(f)
        else:
            self._setPeptideHits(results)
    
    # Internal function that sets peptide hits from a Path object representing a single file.
    def _setPeptideHits(self, resultFile):
        with resultFile.open(mode='r') as tsvin:
            tsvReader = csv.reader(tsvin, delimiter='\t')
            for row in tsvReader:
                protType = determineIDType(row)
                if protType == 'first':
                    continue
                if not isSignificant(row):
                    break
                if protType != 'contaminant' and protType != 'decoy':
                    if not row[PEPTIDE] in self.peps.keys():
                        self.peps[row[PEPTIDE]] = set()
                    for hit in getProteinHitList(row, protType):
                        self.peps[row[PEPTIDE]].add(hit)
    
    # Retrieves a protein from the reference database by its ID.
    def getProt(self, id):
        return self.db[id]
    
    # Returns a list of protein IDs that matched the given peptide sequence
    def getProteinIDsForPeptide(self, sequence):
        return list(self.peps[sequence])
    
    # Returns a list of proteins that matched the given peptide sequence
    def getProteinsForPeptide(self, sequence):
        toReturn = []
        for protID in self.peps[sequence]:
            toReturn.append(self.getProt(protID))
        return toReturn
    
    # Returns a list of taxa that matched the given peptide sequence
    def getTaxaForPeptide(self, sequence):
        toReturn = set()
        for protID in self.peps[sequence]:
            prot = self.getProt(protID)
            for taxa in prot.taxa:
                toReturn.add(taxa)
        return list(toReturn)
    
    # Returns a list of protein names that matched the given peptide sequence
    def getProteinNamesForPeptide(self, sequence):
        toReturn = []
        for prot in self.peps[sequence]:
            toReturn.append(prot.name)
        return toReturn
    
    # Get the list of taxa associated with the protein ID.
    def getTaxaFor(self, id):
        return self.db[id].taxa
    
    # Returns a randomly ordered list of all the taxa with proteins in this reference database
    def getTaxaList(self):
        return list(self.taxa)
    
    # Returns a randomly ordered list of all the genuses with proteins in this reference database
    def getGenusList(self):
        toReturn = set()
        for taxa in self.taxa:
            toReturn.add(getFirstWord(taxa))
        return list(toReturn)
    
    # Creates a new TaxaCounter object with entries for each taxa represented in this ProtRef database
    def generateTaxaCounter(self):
        allTaxa = []
        for x in self.db.values():
            allTaxa.append(x.taxa)
        return TaxaCounter(self)

# Object that stores all the genomes from a ProtRef object and tracks what peptides each covers based on a set of TSV result files.
class GenomeRef:
    def __init__(self, protRef, tsvPath):
        self.db = {}
        for tsvFile in [tsvFile for tsvFile in tsvPath.iterdir() if tsvFile.suffix == '.tsv']:
            with tsvFile.open(mode='r') as infile:
                tsvReader = csv.reader(infile, delimiter='\t')
                for row in tsvReader:
                    if row[Q_VALUE] == 'QValue':
                        continue
                    if float(row[Q_VALUE]) > 0.01:
                        break
                    if determineIDType(row) != 'bacteria':
                        continue
                    hitList = getProteinHitList(row, 'bacteria')
                    for hit in hitList:
                        protRef.getProt(hit).hitPeptides.add(row[PEPTIDE])
                        gens = protRef.getTaxaFor(hit)
                        for genName in gens:
                            if not genName in self.db.keys():
                                self.db[genName] = Genome(genName)
                            self.db[genName].addPeptide(row[PEPTIDE])
    
    # Add a genome to this GenomeRef. Doesn't overwrite if the genome is already in the ref.
    def addGenome(self, name):
        if not name in self.db.keys():
            self.db[name] = Genome(name)
    
    # Returns the genome object with the specified name
    def getGenome(self, name):
        return self.db[name]
    
    # Returns a list of all the genomes in this GenomeRef
    def getGenomes(self):
        toReturn = []
        for genome in self.db.values():
            toReturn.append(genome)
        return toReturn

# Object that stores NGS abundance data for a single sample.
class NGSSample:
    def __init__(self, ID):
        self.data = {}
        self.id = ID
    
    # Associate the percent abundance with the specified taxa for this sample
    def addData(self, taxa, percent):
        self.data[taxa] = percent
    
    # Get the percent abundance for the sample of the specified taxa. If the taxa is not in this sample, returns 0
    def getAbundance(self, taxa):
        if taxa in self.data.keys():
            return self.data[taxa]
        else:
            return 0
    
    # Removes the data for the given taxa from this sample.
    def removeTaxa(self, taxa):
        if taxa in self.data.keys():
            del self.data[taxa]
    
# An object that stores NGS data for different organisms per sample and provides methods to easily access that data.
class NGSRef:
    def __init__(self, csvPath, cutoff):
        self.samples = []
        self.taxa = set()
        with csvPath.open(mode='r') as infile:
            reader = csv.reader(infile)
            headerFlag = True
            for row in reader:
                if headerFlag:
                    headerFlag = False
                    continue
                percent = float(row[3])
                if percent < cutoff:
                    continue
                self.addDataToSample(row[0], row[1], percent)
                self.taxa.add(row[1])
    
    # Note the abundance for the specified taxa in the sample with the specified ID
    def addDataToSample(self, ID, taxa, abundance):
        sample = None
        try:
            sample = self.getSample(ID)
        except:
            sample = NGSSample(ID)
            self.samples.append(sample)
        sample.addData(taxa, abundance)
    
    # Combine the 2 taxa into a single taxa for each sample under 'newName'
    def combineTaxa(self, taxa1, taxa2, newName):
        if not taxa1 in self.taxa:
            raise Exception(f'{taxa1} is not in this reference.')
        if not taxa2 in self.taxa:
            raise Exception(f'{taxa2} is not in this reference.')
        for i in range(len(self.samples)):
            self.samples[i].addData(newName, self.samples[i].getAbundance(taxa1) + self.samples[i].getAbundance(taxa2))
            self.samples[i].removeTaxa(taxa1)
            self.samples[i].removeTaxa(taxa2)
        self.taxa.remove(taxa1)
        self.taxa.remove(taxa2)
        self.taxa.add(newName)
    
    # Get the percent abundance for the taxa in the specified sample
    def getAbundance(self, ID, taxa):
        return self.getSample(ID).getAbundance(taxa)
    
    # Returns a list containing the abundance for the taxa in each sample, in order
    def getAbundanceList(self, taxa):
        toReturn = []
        for sample in self.samples:
            toReturn.append(sample.getAbundance(taxa))
        return toReturn
    
    # Return a list of the names of all taxa present in samples in this NGSRef
    def getTaxaList(self):
        return list(self.taxa)
    
    # Get the sample with the specified ID
    def getSample(self, ID):
        for sample in self.samples:
            if sample.id == ID:
                return sample
        raise Exception(f'No sample with ID {ID}')

# Object used to quickly grab information about GO numbers. Initialized directly from the GO database.
class GoRef:
    def __init__(self, filePath):
        self.db = {}
        rawText = ''
        with filePath.open(mode='r') as infile:
            rawText = infile.read()
        entries = rawText.split('[Term]\n')
        for entry in entries:
            if entry == '':
                continue
            lines = entry.split('\n')
            identifier = lines[0][4:]
            description = lines[1][6:]
            self.db[identifier] = description
    
    # Returns the function for the specified GO number.
    def getFunc(self, go):
        try:
            return self.db[go]
        except:
            if go.find('GO:') != -1:
                return '.'
            else:
                raise Exception(f'{go} could not be found in the database.')
    
    # Returns a list of all the go numbers with the given function. Case insensitive.
    def numbersForFunction(self, funcString):
        funcString = funcString.lower()
        toReturn = []
        for key, value in self.db.items():
            if value.lower() == funcString:
                toReturn.append(key)
        return toReturn

# Object used to track counts of arbitrary stuff for a group of different taxa.
# Uses a hash set as the data backend so accessing counts for each taxa is O(1).
class TaxaCounter:
    def __init__(self, protRef):
        self.prots = protRef
        self.db = {}
        for taxaName in protRef.getTaxaList():
            self.db[taxaName] = 0
    
    # Increment the count for the correct taxa based on protein ID
    def increment(self, id):
        self.db[self.prots.getProt(id).taxa] += 1
    
    # Get the count for the specified taxa
    def getCount(self, taxa):
        return self.db[taxa]
    
    # Returns a list of all the taxa in this counter in alphabetical order
    def getAlphabeticalTaxa(self):
        toReturn = list(self.db.keys())
        toReturn.sort(key=lambda x: x, reverse=False)
        return toReturn
    
    # Returns a list of all the counts for the alphabetical-ordered taxa
    def getAlphabeticalCounts(self):
        tupleHolder = []
        for x in self.db.keys():
            tupleHolder.append((x, self.db[x]))
        tupleHolder.sort(key=lambda x: x[0], reverse=False)
        toReturn = []
        for y in tupleHolder:
            toReturn.append(y[1])
        return toReturn
    
    # Returns the ProtRef database underlying this TaxaCounter
    def getProteinRef(self):
        return self.prots
    
    # Resets the count for each taxa to 0
    def reset(self):
        for x in self.db.keys():
            self.db[x] = 0

# Object used to track counts of proteins for different genuses.
# Identifies a genus by the first word in the taxa name.
# Uses a hash set as the data backend so accessing counts for each taxa is O(1).
class GenusCounter(TaxaCounter):
    def __init__(self, protRef):
        self.prots = protRef
        self.db = {}
        for taxaName in protRef.getTaxaList():
            if taxaName.find('Homo sapiens') != -1:
                continue
            self.db[getFirstWord(taxaName)] = 0
    
    # Increment the count for the correct taxa based on protein ID
    def increment(self, id):
        self.db[getFirstWord(self.prots.getTaxaFor(id))] += 1
    
    # Returns the genus for the specified protein ID
    def getTaxaFor(self, id):
        return getFirstWord(self.prots.getTaxaFor(id))

# Class used for getting the genus, species, and strain name of an organism
class Taxonomy:
    def __init__(self, taxaTxt):
        txt = taxaTxt.split(' ')
        self.genus = txt[0]
        self.species = ''
        self.strain = ''
        if len(txt) > 1:
            self.species = txt[1]
        if len(txt) > 2:
            for i in range(2, len(txt)):
                self.strain += txt[i]
                if i != len(txt) - 1:
                    self.strain += ' '

# Class representing a uniform 2-dimensional array of objects.
class Array2D:
    def __init__(self, data=[], rows=0, cols=0, val=None):
        if data != [] and rows != 0:
            raise Exception('Error: Cannot specify data to create 2D array and attempt to set length.')
        self.rows = []
        for i in range(rows):
            self.rows.append([val] * cols)
        for datum in data:
            self.rows.append(datum)
    
    # Add a new row to the 2D array
    def appendRow(self, newRow):
        if len(self.rows) != 0 and len(self.rows[0]) != len(newRow):
            raise Exception('Row length does not match the current array.')
        self.rows.append(newRow)
    
    # Add a new column to the 2D array
    def appendCol(self, newCol):
        if len(self.rows) == 0:
            for i in range(len(newCol)):
                self.rows.append([])
        if len(newCol) != len(self.rows):
            raise Exception('Column length does not match the current array.')
        for i in range(len(self.rows)):
            self.rows[i].append(newCol[i])
    
    # Returns the row at the specified index in a list
    def getRow(self, index):
        return self.rows[index].copy()
    
    # Returns the column at the specified index in a list
    def getCol(self, index):
        toReturn = []
        for row in self.rows:
            toReturn.append(row[index])
        return toReturn
    
    # Get the value stored in the array in the specified row and column
    def getVal(self, row, col):
        return self.rows[row][col]
    
    # Set the value in the specified row and column of the array
    def setVal(self, row, col, newVal):
        self.rows[row][col] = newVal
    
    # Resets the value of an entire row at the specified index
    def setRow(self, index, newVals):
        if len(newVals) != len(self.rows[index]):
            raise Exception('Error: The number of values in the specified list does not match the number of values in the row to be changed')
        self.rows[index] = newVals
    
    # Resets the value of an entire column at the specified index
    def setCol(self, index, newVals):
        if len(newVals) != len(self.rows):
            raise Exception('Error: Column length does not match the current array.')
        for i in range(len(self.rows)):
            row[i][index] = newVals[i]
    
    # Perform the specified operation on every value in the table
    def operateAll(self, operation):
        for i in range(len(self.rows)):
            for j in range(len(self.rows[i])):
                #print(self.rows[i][j])
                self.rows[i][j] = operation(self.rows[i][j])
                #print(self.rows[i][j])
    
    # Adds each individual value in 'otherArray' to the corresponding value in this 2D array
    def add(self, otherArray):
        if otherArray.length() != self.length():
            raise Exception('Error: 2D arrays must be the same size to add them together')
        for i in range(len(self.rows)):
            for j in range(len(self.rows[i])):
                self.rows[i][j] += otherArray.rows[i][j]
    
    # Returns a tuple where the first value is the length of rows in the 2D array, and the second value is the length of columns
    def length(self):
        if len(self.rows) == 0:
            return (0, 0)
        return (len(self.rows[0]), len(self.rows))

# Class representing a 3-dimensional array.
class Array3D:
    def __init__(self, tables=[]):
        self.tables = []
        for t in tables:
            self.appendTable(t)
    
    # Add a new table to the array
    def appendTable(self, newTable):
        self.tables.append(newTable)
    
    # Get the 2D array from the 3D array at the specified index
    def getTable(self, index):
        return self.tables[index]
    
    # Returns the value at the specified row and index from each 2D array in the 3D array as a list
    def getStack(self, row, col):
        toReturn = []
        for table in self.tables:
            toReturn.append(table.getVal(row, col))
        return toReturn
    
    # Returns a 2D array where the value of each cell is the average of that cell position in all tables in the 3D array
    def getAverage(self):
        toReturn = Array2D(rows=self.tables[0].length()[0], cols=self.tables[0].length()[1], val=0)
        for i in range(self.tables[0].length()[0]):
            for j in range(self.tables[0].length()[1]):
                toReturn.setVal(i, j, sum(self.getStack(i, j)) / self.depth())
        return toReturn
    
    # Returns the number of tables in this 3D array.
    def depth(self):
        return len(self.tables)

# This function takes a row from a tsv result file and returns a list of protein ID's [protein1, protein2, ...]
# If the proteinString is from the first row in the TSV file (proteinString = 'Protein'), the function will return 'Protein'
# 'protType' can be used to filter the type of proteins the function returns. If 'all' is specified, the function returns all proteins hit
def getProteinHitList(row, protType):
    if type(row) == 'str':
        raise Exception('getProteinHitList method takes a row as its first argument.')
    if row[PROTEIN_HITS] == 'Protein':
        return ['Protein']
    proteins = row[PROTEIN_HITS].split(';')
    toReturn = []
    for i in range(0, len(proteins)):
        endIndex = proteins[i].find('(')
        hitID = proteins[i][0:endIndex]
        if protType == 'all' or protTypeIs(hitID, protType):
            toReturn.append(hitID)
    return toReturn

# This function takes a single protein entry from a .fasta file and returns the species name for that entry.
def extractSpeciesName(entry):
    stringStart = entry.find('OS=') + 3
    stringEnd = entry.find(' OX=')
    return entry[stringStart:stringEnd]

# This function takes a single protein entry from a .fasta file and returns the protein identifier for that entry.
def extractProteinID(entry):
    stringEnd = entry.find(' ')
    return entry[0:stringEnd]

# This function takes a single protein entry from a .fasta file and returns the protein identifier for that entry.
def extractProteinName(entry):
    stringStart = entry.find(' ') + 1
    stringEnd = entry.find(' OS=')
    return entry[stringStart:stringEnd]

# Returns the first word in a string
def getFirstWord(text):
    words = text.split(' ')
    if words[0] == '':
        return words[1]
    else:
        return words[0]

# Returns the protein's type depending on what protein set it is in.
# Can return 'human', 'contaminant', 'bacteria', 'fungi', or 'trichomonas'
def determineProteinType(protID):
    if protID in HUMAN_PROTEINS:
        return 'human'
    elif protID in BACTERIAL_PROTEINS:
        return 'bacteria'
    elif protID in CONTAMINANT_PROTEINS:
        return 'contaminant'
    elif protID in FUNGAL_PROTEINS:
        return 'fungi'
    elif protID in TRICHOMONAS_PROTEINS:
        return 'trichomonas'
    else:
        raise Exception(f'Could not identify protein type: "{protID}"')

# Takes a spectral hit and determines the type of protein it represents.
# A spectrum can be matched to proteins from multiple organisms, i.e. both human and bacteria, so takes a hierarchical approach.
# decoy > contaminant > human > fungi > trichomonas > bacteria
# Can return the strings 'human', 'decoy', 'contaminant', 'first' (when the "protein" is actually the first row of the TSV file),
# 'fungi', 'trichomonas', or 'bacteria'
def determineIDType(spectrumRow):
    protType = None
    for prot in getProteinHitList(spectrumRow, 'all'):
        if prot == 'Protein':
            return 'first'
        elif prot.find('XXX') != -1:
            return 'decoy'
        protID = determineProteinType(prot)
        if protID == 'contaminant':
            return 'contaminant'
        elif protID == 'human' and (protType == None or protType == 'bacteria' or protType == 'fungi' or protType == 'trichomonas'):
            protType = 'human'
        elif protID == 'fungi' and (protType == None or protType == 'bacteria' or protType == 'trichomonas'):
            protType = 'fungi'
        elif protID == 'trichomonas' and (protType == None or protType == 'bacteria'):
            protType = 'trichomonas'
        elif protID == 'bacteria' and protType == None:
            protType = 'bacteria'
    if protType == None:
        raise Exception(f'Prottype was None: {str(spectrumRow)}')
    return protType

# Used to verify that the protein with the given ID is of the correct type.
# Generally used to verify that a given protein is the same type as the spectrum has been associated with
def protTypeIs(proteinID, proteinType):
    if proteinType == 'human':
        return proteinID in HUMAN_PROTEINS
    elif proteinType == 'bacteria':
        return proteinID in BACTERIAL_PROTEINS
    elif proteinType == 'contaminant':
        return proteinID in CONTAMINANT_PROTEINS
    elif proteinType == 'decoy':
        return proteinID.find('XXX') != -1
    elif proteinType == 'fungi':
        return proteinID in FUNGAL_PROTEINS
    elif proteinType == 'trichomonas':
        return proteinID in TRICHOMONAS_PROTEINS
    else:
        raise ValueError(f'Illegal protein type "{proteinType}"')

# Takes in a row representing a spectrum hit. Returns True if the Q value is < 0.01, False otherwise.
def isSignificant(row):
    if float(row[Q_VALUE]) < 0.01:
        return True
    return False

# Determines if the spectrum hit on the specified row matches stringent criteria:
# Length of the identified peptide must be at least 10 amino acids
# The peptide only matched one taxa in the database
def isStringent(spectrumRow, protRef, genus=False):
    if type(spectrumRow) == 'str':
        raise Exception('isStringent method takes a row as its first argument.')
    if len(spectrumRow[PEPTIDE].replace('+15.995', '')) < 10:
        return False
    hits = getProteinHitList(spectrumRow, 'all')
    if len(hits) == 1:
        return True
    firstTaxaList = protRef.getTaxaFor(hits[0])
    firstTaxaHit = firstTaxaList[0] if genus == False else getFirstWord(firstTaxaList[0])
    if len(firstTaxaList) > 1:
        for i in range(1, len(firstTaxaList)):
            newTaxa = firstTaxaList[i] if genus == False else getFirstWord(firstTaxaList[i])
            if newTaxa != firstTaxaHit:
                return false
    for i in range(1, len(hits)):
        nextTaxaList = protRef.getTaxaFor(hits[i])
        for i in range(len(nextTaxaList)):
            newTaxa = nextTaxaList[i] if genus == False else getFirstWord(nextTaxaList[i])
            if newTaxa != firstTaxaHit:
                return False
    return True

# This function is used to parse through a set of TSV files, extracting data from each one.
# The function takes in a string representing the relative path to the folder containing the TSV files.
# The function also needs a method which takes a csv.reader() object, parses through each row of the file, and returns some value.
# This function then adds the returned value to a generic list, which it returns at the end of method execution.
def parseTSVs(directoryPath, parseMethod):
    toReturn = []
    for file in directoryPath.iterdir():
        if file.suffix != '.tsv':
            continue
        with file.open(mode='r') as tsvin:
            tsvReader = csv.reader(tsvin, delimiter='\t')
            toReturn.append(parseMethod(tsvReader))
    return toReturn

# This function is used to parse through a set of TSV files, extracting data from each one.
# The function takes in a string representing the relative path to the folder containing the TSV files.
# The function also needs a method which takes a csv.reader() object as well as a list of arguments, parses through each row of the file, and returns some value.
# This function then adds the returned value to a generic list, which it returns at the end of method execution.
def parseTSVWithArgs(directoryPath, parseMethod, parseArgs):
    toReturn = []
    for file in directoryPath.iterdir():
        if file.suffix != '.tsv':
            continue
        with file.open(mode='r') as tsvin:
            tsvReader = csv.reader(tsvin, delimiter='\t')
            toReturn.append(parseMethod(tsvReader, parseArgs))
    return toReturn

# Returns True if 'item' is inside of 'searchList'.
# Note: 'item' must have an equals() method
def itemInList(item, searchList):
    for x in searchList:
        if item.equals(x):
            return True
    return False

# Returns a set of peptides that pass protein filtering criteria
# Two unique peptides that match to the same protein, peptides with two spectral hits, or spectral probability < 1e-15
# typeOfProt must be either 'human', 'bacteria', 'fungi', or 'trichomonas'
def getFilteredPeptides(results, typeOfProt):
    proteins = {} # key=protein id, value=set of peptide sequences
    pepCounts = {} # key=peptide sequence, value=number of times identified across all samples
    pepProbabilities = {} # key=peptide sequence, value=lowest spectral probability
    if not isinstance(results, list):
        results = [x for x in results.iterdir() if x.suffix == '.tsv']
        results.sort(key=lambda x: x.name)
    for res in results:
        with res.open(mode='r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for row in reader:
                if row[PROTEIN_HITS] == 'Protein':
                    continue
                if not isSignificant(row):
                    break
                protType = determineIDType(row)
                if protType == typeOfProt:
                    hits = getProteinHitList(row, protType)
                    for hit in hits:
                        if not hit in proteins.keys():
                            proteins[hit] = set()
                        proteins[hit].add(row[PEPTIDE])
                    if not row[PEPTIDE] in pepCounts.keys():
                        pepCounts[row[PEPTIDE]] = 0
                    pepCounts[row[PEPTIDE]] += 1
                    if not row[PEPTIDE] in pepProbabilities.keys():
                        pepProbabilities[row[PEPTIDE]] = 100000000000000000000
                    prob = float(row[SPEC_PROBABILITY])
                    if prob < pepProbabilities[row[PEPTIDE]]:
                        pepProbabilities[row[PEPTIDE]] = prob
    passingPeps = set()
    for peps in proteins.values():
        if len(peps) > 1:
            for pep in peps:
                passingPeps.add(pep)
    for pep, probability in pepProbabilities.items():
        if probability < 1e-15:
            passingPeps.add(pep)
    return passingPeps


# Returns a list of lists where each sub-list contains all the peptides from 'allowedPeps' for the result file, in order, with repeats
def getPeptidesByResult(results, allowedPeps):
    toReturn = []
    for res in results:
        toAppend = []
        with res.open(mode='r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for row in reader:
                protType = determineIDType(row)
                if protType == 'first':
                    continue
                if not isSignificant(row):
                    break
                if row[PEPTIDE] in allowedPeps:
                    toAppend.append(row[PEPTIDE])
        toReturn.append(toAppend)
    return toReturn

# Writes head, contents of listToWrite to the outputFile
def writeListToFile(head, listToWrite, outputFile):
    # Write the results to the output file
    with outputFile.open(mode='w') as outfile:
        outfile.write(head)
        outfile.write(''.join(listToWrite))
    print('Results written to file')

# Returns the header string for a func file
def getFuncHeader():
    return 'sample\ttaxa\toriginal peptide\tprotein\tpeptide queried\ttotal_protein_count\tgo_term\n'

### GRAPHING FUNCTIONS ###
# Creates a normalized stacked bar chart for the specified data.
# 'samples' are the individual bar groups - i.e. the labels that go on the x axis
# 'categories' are the groups that build a single bar - i.e. the values that will go in the legend
# 'data' is a 2D array structured as: [[sample1A, sample1B], [sample2A, sample2B], [sample3A, sample3B]]
# If 'sort' is True, the function sorts the categories by alphabetical order before graphing them.
def generateNormalizedStackedBar(samples, categories, data, title='', ylabel='Percent Detected Spectra', xLabTotals=False, sort=False, colors=COLORS, savepath=None, legendOrder=None):
    toGraph = data.copy()
    tempCategories = categories.copy()
    if sort:
        tempData = Array2D(data=toGraph)
        tupleHolder = []
        for i in range(len(categories)):
            tupleHolder.append((categories[i], tempData.getCol(i)))
        tupleHolder.sort(key=lambda x: x[0], reverse=False)
        tempSorted = Array2D()
        tempCategories = []
        for tup in tupleHolder:
            tempSorted.appendCol(tup[1])
            tempCategories.append(tup[0])
        toGraph = tempSorted.rows
    totals = []
    for row in toGraph:
        totals.append(sum(row))
    normalized = []
    for i in range(len(toGraph[0])):
        normRow = []
        for j in range(len(toGraph)):
            if totals[j] == 0:
                normRow.append(0)
            else:
                normRow.append(toGraph[j][i] / totals[j] * 100)
        normalized.append(normRow)
    bottoms = [0] * len(normalized[0])
    for i in range(len(normalized)):
        plt.bar(list(range(len(normalized[0]))), normalized[i], bottom=bottoms, color=colors[i], edgecolor='white', label=tempCategories[i])
        for j in range(len(normalized[i])):
            bottoms[j] += normalized[i][j]
    if xLabTotals:
        plt.xticks(list(range(len(normalized[0]))), totals)
        plt.xlabel('Total Spectra in Sample')
    else:
        plt.xticks(list(range(len(normalized[0]))), samples)
    if not title == '':
        plt.title(title)
    plt.ylabel(ylabel)
    if not legendOrder == None:
        order = legendOrder
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper left', bbox_to_anchor=(1,1), ncol=math.ceil(len(tempCategories)/30))
    else:
        plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=math.ceil(len(tempCategories)/30))
    if not savepath == None:
        plt.savefig(savepath, bbox_inches='tight', dpi=300)

# Creates a grouped bar graph where 'xCategories' are the categories along the x-axis the data falls into, 'groupsLabel' are what to call the groups in the legend,
# 'groups' is a list of strings defining the different groups, and data is a list of lists where data[0] is the data for groups[0], etcl.
# Returns the resultant bar graph object.
def generateGroupedBar(xCategories, groupsLabel, groups, data, xlabel='', ylabel='', title='', size=(6, 2)):
    dataList = []
    for i in range(len(data)):
        dataList.extend(data[i])
    catList = []
    groupList = []
    for i in range(len(groups)):
        catList += xCategories
        for j in range(len(data[0])):
            groupList.append(groups[i])
    dataFrame = pd.DataFrame({'data':dataList, 'xCategories':catList, groupsLabel:groupList})
    graph = sns.catplot(x='xCategories', y='data', hue=groupsLabel, data=dataFrame, height=size[0], aspect=size[1], kind="bar")
    graph.despine(left=True)
    graph.set_ylabels(ylabel)
    graph.set_xlabels(xlabel)
    graph.fig.suptitle(title)
    return graph

# Convenience function for setting various parameters in pyplot.
def setPlt(size=(0, 0), xlim=(0, 0), ylim=(0, 0), xlabel='', ylabel='', title='', xticksrot=0):
    if size != (0, 0):
        plt.figure(figsize=size)
    if xlim != (0, 0):
        plt.xlim(xlim[0], xlim[1])
    if ylim != (0, 0):
        plt.ylim(ylim[0], ylim[1])
    if xlabel != '':
        plt.xlabel(xlabel)
    if ylabel != '':
        plt.ylabel(ylabel)
    if title != '':
        plt.title(title)
    if xticksrot != 0:
        plt.xticks(rotation=xticksrot, ha='right')
    plt.tight_layout()

# Convenience function that saves a 300DPI PNG of the specified figure to the specified filePath
def saveFig(figure, filePath):
    figure.get_figure().savefig(filePath, bbox_inches='tight', dpi=300)

# Takes in a file of amino acid sequences and outputs a new file with the tryptic peptides resulting from those sequences.
def prot2pept(sequenceFile, peptideFile):
    trypticPeptides = set()
    with sequenceFile.open(mode='r') as infile:
        for sequence in infile:
            lastCleavage = 0
            for i in range(len(sequence)):
                if sequence[i] == '\n':
                    trypticPeptides.add(sequence[lastCleavage:i])
                    newPeptide = False
                elif sequence[i] == 'R' or sequence[i] == 'K':
                    if sequence[i + 1] == 'P':
                        continue
                    trypticPeptides.add(sequence[lastCleavage:i + 1])
                    lastCleavage = i + 1
    with peptideFile.open(mode='w', newline='') as outfile:
        for pep in trypticPeptides:
            if len(pep) > 4 and len(pep) < 50:
                outfile.write(f'{pep}\n')

# Runs the unipept pept2go function on all the peptides in the 'peptideFile', then writes the results to 'outputFile'.
def pept2go(peptideFile, outputFile):
    header = ''
    toWrite = ''
    unipept = shutil.which('unipept')
    count = 0
    for pep in peptideFile.read_text().split('\n'):
        if pep == '':
            continue
        command = [unipept, 'pept2go', pep, '--all']
        completedProcess = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        if completedProcess.stdout == b'':
            toWrite = f'{toWrite}{pep},-,-,-,-\n'
            continue
        processTxt = completedProcess.stdout.decode('UTF-8').split('\r\n')
        if header == '':
            header = processTxt[0] + '\n'
        toWrite = f'{toWrite}{processTxt[1]}\n'
        count += 1
        if count % 50 == 0:
            print(f'{str(count)} peptides annotated.')
    with outputFile.open(mode='w', newline='') as outfile:
        outfile.write(header)
        outfile.write(toWrite)

# Takes in a string representing an amino acid sequence that may or may not have newline or space characters in it and returns a list of tryptic peptides that would result
def getTrypticPeptides(aaString):
    aaString = aaString.replace('\n', '')
    aaString = aaString.replace(' ', '')
    tryptics = []
    lastCleavage = 0
    for i in range(len(aaString)):
        if i == len(aaString) - 1:
            tryptics.append(aaString[lastCleavage:i + 1])
        elif aaString[i] == 'R' or aaString[i] == 'K':
            if aaString[i + 1] == 'P':
                continue
            tryptics.append(aaString[lastCleavage:i + 1])
            lastCleavage = i + 1
    return tryptics

# Returns the longest tryptic peptide that is <= 30 amino acids long. If there are no tryptics peptides <= 30 aa's returns the smallest.
def getBestTrypticPeptide(aaString):
    tryptics = getTrypticPeptides(aaString)
    toQuery = []
    if len(tryptics) > 2:
        toQuery = tryptics[1:len(tryptics) - 2]
    else:
        return tryptics[0]
    toQuery.sort(key=lambda x: len(x), reverse=True)
    for pep in toQuery:
        if len(pep) <= 30:
            return pep
    return tryptics[len(tryptics) - 1]

# Get functions for all peptides in a list of peptide lists.
# ProtRefs can be a list of protrefs to work with Tailored database results, or a single ProtRef for Community database results
# pepType can be 'human' or 'bacteria'
# Returns a large list of data that is ready to be written out to file.
def peps2func(listPeptidesList, protRefs, pepType, waitTime=3):
    # Query unipept for function information
    data = {}
    toWrite = []
    unipept = shutil.which('unipept')
    count = 0
    totalPeps = 0
    _iter = -1
    for i in range(len(listPeptidesList)):
        totalPeps += len(listPeptidesList[i])
    for pepList in listPeptidesList:
        _iter += 1
        ref = None
        if type(protRefs) == list:
            ref = protRefs[_iter]
        else:
            ref = protRefs
        for pep in pepList:
            if not pep in data.keys():
                taxaString = 'ERROR'
                if pepType == 'human':
                    taxaString = 'Homo sapiens'
                elif pepType == 'bacteria':
                    taxaList = ref.getTaxaForPeptide(pep)
                    for taxa in taxaList:
                        taxaString += taxa + ';'
                if taxaString == 'ERROR':
                    raise Exception('Invalid value for pepType')
                matchingProtein = ref.getProteinIDsForPeptide(pep)[0]
                trypticPeps = getTrypticPeptides(pep.replace('+15.995', ''))
                toQuery = ''
                if len(trypticPeps) == 1:
                    toQuery = trypticPeps[0]
                elif len(trypticPeps) == 2:
                    lastChar = trypticPeps[1][len(trypticPeps[1]) - 1]
                    if lastChar == 'K' or lastChar == 'R':
                        toQuery = trypticPeps[1]
                    else:
                        toQuery = trypticPeps[0] if len(trypticPeps[0]) > len(trypticPeps[1]) else trypticPeps[1]
                else:
                    lastPep = trypticPeps[len(trypticPeps) - 1]
                    lastChar = lastPep[len(lastPep) - 1]
                    searchIndex = len(trypticPeps) if lastChar == 'K' or lastChar == 'R' else len(trypticPeps) - 1
                    fullyTryptic = trypticPeps[1:searchIndex]
                    for trypPep in fullyTryptic:
                        if len(trypPep) > len(toQuery):
                            toQuery = trypPep
                command = [unipept, 'pept2go', toQuery, '--all']
                toPrint = ''
                try:
                    completedProcess = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True, check=True, timeout=10)
                    if completedProcess.stdout == b'':
                        raise Exception('Go to empty')
                    processTxt = completedProcess.stdout.decode('UTF-8').split('\r\n')
                    queryOutput = processTxt[1].split(',')
                    if len(queryOutput) < 4:
                        raise Exception('Go to empty')
                    goList = queryOutput[:3]
                    goData = f'{goList[0]}\t{goList[1]}\t{goList[2]}'
                    comm = ','
                    sep = '\t'
                    data[pep] = f'{taxaString}\t{pep}\t{matchingProtein}\t{goData.replace(comm, sep)}\n'
                    toPrint = f'Got response for {toQuery} {str(count + 1)}/{str(totalPeps)}'
                except:
                    data[pep] = f'{taxaString}\t{pep}\t{matchingProtein}\t{toQuery}\t-\t-\n'
                    toPrint = f'Had to except {toQuery} {str(count + 1)}/{str(totalPeps)}'
            else:
                toPrint = f'{toQuery} already queried. {str(count + 1)}/{str(totalPeps)}'
            count += 1
            if count % 1000 == 0:
                print(toPrint)
            toWrite.append(f'{SAMPLE_NAMES[_iter]}\t{data[pep]}')
    return toWrite

# This function condenses the results from a large, broken-up database search into a single, stitched-together TSVFile object.
# Matches result files together by the first 4 characters in their file names.
# Decides which result is best based on which has the highest MSGFscore.
# Returns the results as a dictionary of TSVFiles
def condenseHugeDBResults(resultsDir):
    fileSorter = {} # key=file id, value=list of files with the file id
    for res in [x for x in resultsDir.iterdir() if x.suffix == '.tsv']:
        fileID = res.name[0:4]
        if not fileID in fileSorter.keys():
            fileSorter[fileID] = []
        fileSorter[fileID].append(res)
    print('Done sorting files')
    toReturn = {} # key=file id, value=TSVFile of stitched-together result
    for fileID, resultFiles in fileSorter.items():
        masterFile = TSVFile(resultFiles[0])
        for i in range(1, len(resultFiles)):
            nextFile = TSVFile(resultFiles[i])
            for specID, row in nextFile.rows.items():
                if not specID in masterFile.rows.keys():
                    masterFile.setRow(specID, row)
                currentScore = masterFile.getMSGFScore(specID)
                newScore = nextFile.getMSGFScore(specID)
                if newScore > currentScore:
                    masterFile.setRow(specID, row)
                elif newScore == currentScore:
                    rowToUpdate = masterFile.getRow(specID)
                    rowToUpdate[PROTEIN_HITS] = rowToUpdate[PROTEIN_HITS] + ';' + row[PROTEIN_HITS]
        print(f'Done condensing file {fileID}')
        toReturn[fileID] = masterFile
    return toReturn

# Returns the files in the given directory of the given suffix in a list, ordered by their name
def getOrderedFiles(directory, suffix):
    toReturn = [x for x in directory.iterdir() if x.suffix == suffix]
    toReturn.sort(key=lambda x: x.name)
    return toReturn

# Reads in all protein sequences in 'directory', and outputs the new file to 'dbName'
# Collapses identical IDs by default, or by identical amino acid sequences if collapseIDs == False
def buildDatabase(directory, dbName, collapseIDs=True):
    toWrite = {} # key=protein ID, value=protein object
    for f in [x for x in directory.iterdir() if x.suffix == '.fasta']:
        with f.open(mode='r') as database:
            entry = []
            completeEntries = []
            for line in database:
                if line[0] == '>' and entry != []:
                    entry[0] = entry[0].replace('>', '')
                    completeEntries.append(''.join(entry))
                    entry = [line]
                else:
                    entry.append(line)
            entry[0] = entry[0].replace('>', '')
            completeEntries.append(''.join(entry))
            for e in completeEntries:
                newProt = Protein(e)
                identifier = newProt.id if collapseIDs else newProt.sequence
                if identifier in toWrite.keys():
                    for t in newProt.taxa:
                        toWrite[identifier].addTaxa(t)
                else:
                    toWrite[identifier] = newProt
    with directory.joinpath(dbName).open(mode='w', newline='') as outfile:
        for prot in toWrite.values():
            outfile.write(prot.getEntry().rstrip() + '\n')
    print(f'{str(len(toWrite))} sequences written.')

# Returns a set of every protein ID that was hit in each TSV file in the directory, ignoring significance
def getHitsInResults(resultsDir):
    protsToInclude = set()
    for res in [x for x in resultsDir.iterdir() if x.suffix == '.tsv']:
        with res.open(mode='r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for row in reader:
                protType = determineIDType(row)
                if protType == 'first':
                    continue
                hitProts = getProteinHitList(row, 'all')
                for hit in hitProts:
                    if hit.find('XXX_') == -1: # Ignore decoy hits
                        protsToInclude.add(hit)
    return protsToInclude

# Reads through every protein in database files in dbFilesDir and, if the ID matches one in prots, pulls it out to write
# to the output database. Collapses proteins with the same ID. Names output database files [baseName]_1.fasta, etc.
def refineHugeDatabase(protsToInclude, dbFilesDir, outputPath, baseName):
    prots = protsToInclude.copy()
    toWrite = {} # key=protein ID, value=protein object
    for f in [x for x in dbFilesDir.iterdir() if x.suffix == '.fasta']:
        with f.open(mode='r') as database:
            entry = []
            for line in database:
                if line[0] == '>' and entry != []:
                    entry[0] = entry[0].replace('>', '')
                    firstLine = entry[0]
                    newID = firstLine[:firstLine.find(' ')]
                    if newID in prots:
                        prots.remove(newID)
                        needEdit = entry[0].find('[') != -1 or entry[0].find('OX=') == -1
                        if needEdit:
                            entry[0] = entry[0].replace('[[', 'OS=')
                            entry[0] = entry[0].replace('[', 'OS=')
                            entry[0] = entry[0].replace(']', '')
                            if entry[0].find('OX=') == -1:
                                entry[0] = entry[0].replace('\n', ' OX=x\n')
                        newProt = Protein(''.join(entry))
                        if newProt.id in toWrite.keys():
                            for t in newProt.taxa:
                                toWrite[newProt.id].addTaxa(t)
                        else:
                            toWrite[newProt.id] = newProt
                    entry = [line]
                else:
                    entry.append(line)
            entry[0] = entry[0].replace('>', '')
            firstLine = entry[0]
            newID = firstLine[:firstLine.find(' ')]
            if newID in prots:
                prots.remove(newID)
                needEdit = entry[0].find('[') != -1 or entry[0].find('OX=') == -1
                if needEdit:
                    entry[0] = entry[0].replace('[[', 'OS=')
                    entry[0] = entry[0].replace('[', 'OS=')
                    entry[0] = entry[0].replace(']', '')
                    if entry[0].find('OX=') == -1:
                        entry[0] = entry[0].replace('\n', ' OX=x\n')
                newProt = Protein(''.join(entry))
                if newProt.id in toWrite.keys():
                    for t in newProt.taxa:
                        toWrite[newProt.id].addTaxa(t)
                else:
                    toWrite[newProt.id] = newProt
        print(f'Done reading {f.name}')
    if len(prots) > 0:
        print(f'Warning: {str(len(prots))} sequences could not be found in the supplied databases.')
    dbCount = 1
    onProt = 0
    protIDs = list(toWrite.keys())
    while onProt < len(protIDs):
        dbPath = outputPath.joinpath(f'{baseName}_{str(dbCount)}.fasta')
        with dbPath.open(mode='a', newline='') as outfile:
            for i in range(onProt, len(protIDs)):
                outfile.write(toWrite[protIDs[i]].getEntry().rstrip() + '\n')
                onProt += 1
                if dbPath.stat().st_size > 300000000:
                    dbCount += 1
                    dbPath = outputPath.joinpath(f'{baseName}_{str(dbCount)}.fasta')
                    break
    print(f'{str(len(toWrite))} sequences written.')

# Takes in result files that may have multiple rows with the same specID because of spectral ambiguity and makes a new
# TSV file for each, collapsing the repeat rows into a single row, then writes the new file out to outputDir. Accomplishes this by
# taking the peptide call and other data from the first row, and appending all protein hits from subsequent repeat rows.
def collapseRepeatRows(results, outputDir):
    for res in results:
        head = None
        rows = {} # key=specID, value=row
        with res.open(mode='r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for row in reader:
                if head == None:
                    head = row
                    continue
                if row[SPEC_ID] in rows.keys():
                    oldRow = rows[row[SPEC_ID]]
                    oldRow[PROTEIN_HITS] = oldRow[PROTEIN_HITS] + ';' + row[PROTEIN_HITS]
                    rows[row[SPEC_ID]] = oldRow
                else:
                    rows[row[SPEC_ID]] = row
        newFile = TSVFile()
        newFile.header = head
        for specid, row in rows.items():
            newFile.setRow(specid, row)
        newFile.writeToFile(outputDir.joinpath(res.name))
        
                    
    











