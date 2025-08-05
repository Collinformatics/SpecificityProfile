import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
import gzip
import os.path
import warnings



class WebApp:
    def __init__(self):
        # Params: Dataset
        self.enzymeName = ''
        self.seqLength = False
        self.subsExp = {}
        self.subsBg = {}

        # Params: Files
        self.fileBg = []
        self.seqExp = None
        self.fileBg = []
        self.seqBg = None
        self.pathData = 'data'
        self.pathSeqs = os.path.join(self.pathData, 'sequences')
        self.pathFigs = os.path.join(self.pathData, 'figures')
        self.pathLog = os.path.join(self.pathData, 'log.txt')

        # Params: Process dna
        self.seq5Prime = False
        self.seq3Prime = False
        self.minPhred = False

        # Params: Filter Dataset
        self.filterPos = False
        self.fixAA = {}
        self.exclAA = {}

        # Params: Figures
        self.datasetTag = None
        self.datasetTagMotif = None
        self.title = ''
        self.titleCombined = ''
        self.titleReleased = ''
        self.titleWeblogo = ''
        self.titleWeblogoCombined = ''
        self.titleWeblogoReleased = ''
        self.titleWords = ''
        self.titleWordsCombined = ''
        # self.figEMSquares = figEMSquares
        # if figEMSquares:
        #     self.figSizeEM = (5, 8)  # (width, height)
        # else:
        self.figSizeEM = (9.5, 8)
        self.figSize = (9.5, 8)
        self.figSizeMini = (self.figSize[0], 6)
        self.residueLabelType = 2  # 0 = full AA name, 1 = 3-letter code, 2 = 1 letter
        self.labelSizeTitle = 18
        self.labelSizeAxis = 16
        self.labelSizeTicks = 13
        self.lineThickness = 1.5
        self.tickLength = 4
        self.colorsAA = self.residueColors()
        self.residues = self.colorsAA.keys()

        # # Params:
        # self. = False
        # self. = False
        # self. = False
        #
        # # Params:
        # self. = False
        # self. = False
        # self. = False


        # Verify directory paths
        if self.pathData is not None:
            if not os.path.exists(self.pathData):
                os.makedirs(self.pathData, exist_ok=True)
        if self.pathSeqs is not None:
            if not os.path.exists(self.pathSeqs):
                os.makedirs(self.pathSeqs, exist_ok=True)
        if self.pathFigs is not None:
            if not os.path.exists(self.pathFigs):
                os.makedirs(self.pathFigs, exist_ok=True)

    @staticmethod
    def residueColors():
        color = ['darkgreen', 'firebrick', 'deepskyblue', 'pink', 'navy', 'black', 'gold']
        # Aliphatic, Acidic, Basic, Hydroxyl, Amide, Aromatic, Sulfur

        return {
            'A': color[0],
            'R': color[2],
            'N': color[4],
            'D': color[1],
            'C': color[6],
            'E': color[1],
            'Q': color[4],
            'G': color[0],
            'H': color[2],
            'I': color[0],
            'L': color[0],
            'K': color[2],
            'M': color[6],
            'F': color[5],
            'P': color[0],
            'S': color[3],
            'T': color[3],
            'W': color[5],
            'Y': color[5],
            'V': color[0]
        }

    def pressButton(self, message):
        print(f'Received data: {message}')

        return {'key': 'Returned data'}

    def filterAA(self):
        print('Processing Substrates')

        return {'AA': 'VEHTVALKQNR'}

    def filterMotif(self):
        print('Filtering Substrates')

        return {'Motif': 'TVALK'}



    def getDatasetTag(self):
        tagFix = 'Fix '
        tagExcl = 'Excl '
        fixNPos = len(self.fixAA)

        if self.exclAA:
            print('Exclude AA:')
            for index, (pos, AA) in enumerate(self.exclAA.items()):
                if len(AA) > 1:
                    tag = f'[{','.join(AA)}]@{pos.replace('fix', '')}'
                else:
                    tag = f'{AA}@{pos.replace('fix', '')}'

                if fixNPos > 1 and index != fixNPos - 1:
                    tagExcl += f'{tag}_'
                else:
                    tagExcl += tag
                print(f'Tag: {tagExcl}')
            self.datasetTag = tagExcl

        if self.fixAA:
            print('Fix AA')
            for index, (pos, AA) in enumerate(self.fixAA.items()):
                if len(AA) > 1:
                    tag = f'[{','.join(AA)}]@{pos.replace('fix', '')}'
                else:
                    tag = f'{AA}@{pos.replace('fix', '')}'

                if fixNPos > 1 and index != fixNPos - 1:
                    tagFix += f'{tag}_'
                else:
                    tagFix += tag
                print(f'Tag: {tagFix}')
            self.datasetTag = tagFix

        self.log(f'Dataset Tag: {self.datasetTag}\n\n')



    def getFilter(self, form):
        self.filterPos = form.get('filterPos', [])
        self.fixAA = {}
        self.exclAA = {}

        # Fix AA
        if any(key.startswith('fix') for key in form.keys()):
            for key, value in form.items():
                if 'fix' in key:
                    self.fixAA[key] = value

            self.fixAA = dict(sorted(self.fixAA.items()))
            print(f'Fixing AA:')
            for key, value in self.fixAA.items():
                print(f'     {key}: {value}')
            print()

        # Exclude AA
        if any(key.startswith('excl') for key in form.keys()):
            for key, value in form.items():
                if 'fix' in key:
                    self.exclAA[key] = value

            self.exclAA = dict(sorted(self.exclAA.items()))
            print(f'Excluding AA:')
            for key, value in self.exclAA.items():
                print(f'     {key}: {value}')
            print()
            print('\nStopping job at: getFilter()\n')
            sys.exit()

        self.getDatasetTag()



    def log(self, txt=None):
        if txt is None:
            with open(self.pathLog, 'w'):
                pass
        else:
            with open(self.pathLog, 'a') as log:
                log.write(f'{txt}\n')



    def evalDNA(self, form):
        self.log()  # Clear the log

        # Process job parameters
        self.enzymeName = form['enzymeName']
        self.fileExp = form['fileExp']
        self.fileBg = form['fileBg']
        self.seq5Prime = form['seq5Prime']
        self.seq3Prime = form['seq3Prime']
        self.seqLength = int(form['seqLength'])
        self.minPhred = int(form['minPhred'])

        # Log job params
        self.log(f'================================= Process DNA '
                 f'==================================')
        self.log(f'5\' Sequence: {self.seq5Prime}\n'
                 f'3\' Sequence: {self.seq3Prime}\n'
                 f'Sequence Length: {self.seqLength}\n'
                 f'Min Phred Score: {self.minPhred}')

        # Evaluate input form
        self.getFilter(form)

        # Params: Tmp
        self.fileExp = 'data/variantsExp.fastq.gz' # Update when using files
        self.fileBg = 'data/variantsBg.fasta'
        # if type(file) == 'fastq':
        #     useQS = True

        # Load the data
        if self.fileExp:
            self.loadDNA(path=self.fileExp, datasetType='Experimental')
        if self.fileBg:
            self.loadDNA(path=self.fileBg, datasetType='Background')

        return {'seq': 'GTGGAACATACCGTGGCGCTGAAACAGAACCGC'}



    def loadDNA(self, path, datasetType):
        # Open the file
        openFn = gzip.open if path.endswith('.gz') else open # Define open function
        with openFn(path, 'rt') as file: # 'rt' = read text mode
            if '.fastq' in path or '.fq' in path:
                data = SeqIO.parse(file, 'fastq')
                warnings.simplefilter('ignore', BiopythonWarning)
            elif '.fasta' in path or '.fa' in path:
                data = SeqIO.parse(file, 'fasta')
                warnings.simplefilter('ignore', BiopythonWarning)
            else:
                self.log(f'ERROR: Unrecognized file\n     {path}\n\n')

            # Translate the dna
            self.translate(data, datasetType, True)



    def translate(self, data, datasetType, forwardRead):
        substrates = {}
        useQS = False
        totalSubsExtracted = 0
        totalSeqsDNA = 0
        printN = 10

        def extractionEfficiency():
            perExtracted = (totalSubsExtracted / totalSeqsDNA) * 100
            self.log(f'Extracted Sequences: {totalSubsExtracted:,}\n'
                     f'Evaluated DNA Sequences: {totalSeqsDNA:,}\n'
                     f'     Extraction Efficiency: {perExtracted} %')

        if forwardRead:
            self.log('Translating dna: Forward Read')
        else:
            self.log('Translating dna: Reverse Read')
        self.log(f'     Dataset: {datasetType}\n')

        # Inspect the file
        for datapoint in data:
            print(f'{datapoint}\n')
            if 'phred_quality' in datapoint.letter_annotations:
                useQS = True
            break

        # Translate DNA
        if useQS:
            for index, datapoint in enumerate(data):
                if totalSubsExtracted >= printN:
                    break

                # Process datapoint
                totalSeqsDNA += 1
                dna = str(datapoint.seq)
                self.log(f'DNA Seq: {dna}')

                # Inspect full dna seq
                if self.seq5Prime in dna and self.seq3Prime in dna:
                    qs = datapoint.letter_annotations['phred_quality']

                    # Find: Substrate indices
                    start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                    end = dna.find(self.seq3Prime)

                    # Extract substrate dna seq
                    substrateDNA = dna[start:end].strip()
                    self.log(f'    Sub: {substrateDNA}')
                    if len(substrateDNA) == self.seqLength * 3:
                        # Express substrate
                        substrate = str(Seq.translate(substrateDNA))
                        self.log(f'    Sub: {substrate}')

                        # Inspect substrate seq: PRINT ONLY
                        if 'X' not in substrate and '*' not in substrate:
                            qs = qs[start:end]
                            self.log(f'     QS: {qs}')
                            if all(score >= self.minPhred for score in qs):
                                self.log(f'Keep Substrate\n')
                                if substrate in substrates.keys():
                                    substrates[substrate] += 1
                                else:
                                    substrates[substrate] = 1
                                totalSubsExtracted += 1
                            else:
                                self.log('')
        else:
            print('Code Me!')
            sys.exit()

        # Evaluate dataquality
        extractionEfficiency()

        sys.exit()
