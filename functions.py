import gzip
import os.path



class WebApp:
    def __init__(self):
        # Params: Dataset
        self.enzymeName = ''
        self.seqLength = False

        # Params: Files
        self.fileBg = []
        self.seqExp = None
        self.fileBg = []
        self.seqBg = None
        self.pathData = 'data'
        self.pathSeqs = os.path.join(self.pathData, 'sequences')
        self.pathFigs = os.path.join(self.pathData, 'figures')
        self.path = ''

        # Params: Process DNA
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



    def evalDNA(self, data):
        print(f'Received Data:')
        for key, value in data.items():
            print(f'     {key}: {value}')
        print('\n')

        # Process job parameters
        self.enzymeName = data['enzymeName']
        self.fileExp = data['fileExp']
        self.fileBg = data['fileBg']
        self.seq5Prime = data['seq5Prime']
        self.seq3Prime = data['seq5Prime']
        self.seqLength = data['seqLength']
        self.filterPos = data.get('filterPos', [])
        self.minPhred = data['minPhred']

        for key, value in data.items():
            if 'fix' in key:
                self.fixAA[key] = value
        self.fixAA = dict(sorted(self.fixAA.items()))

        if 'fix' in data.keys():
            print(f'Fixing AA:')
            for key, value in self.fixAA.items():
                print(f'     {key}: {value}')
            print()


        # Params: Tmp
        # self.fileExp = 'data/variantsExp.fastq.gz'
        # self.fileBg = 'data/variantsBg.fasta.zip'

        # Load the data
        if self.fileExp:
            self.loadDNA(path=self.fileExp)
        if self.fileBg:
            self.loadDNA(path=self.fileBg, loadExp=False)

        return {'seq': 'GTGGAACATACCGTGGCGCTGAAACAGAACCGC'}



    def loadDNA(self, path, loadExp=True):
        print(f'================================= Loading Data '
              f'=================================')
        print(f'File Path:\n      {path}\n')

        # Define open function
        openFn = gzip.open if path.endswith('.gz') else open
        with openFn(path, 'rt') as file:  # 'rt' = read text mode
            for index, line in enumerate(file):
                print(f'      {line.strip()}')
                if index == 10:
                    break
        print('\n')

        # if '.fastq' in path:
        #     print(f'Loading: fastq')
        #     with open(path) as dna:
        #         for line in dna:
        #             print(f'      {line}')
        #     print('\n')
        # elif '.fasta' in path:
        #     print(f'Loading: fasta')
        #     with open(path) as dna:
        #         print(dna)
        #     print('\n')

