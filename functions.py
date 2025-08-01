import os.path


class WebApp:
    def __init__(self):
        # Params: Enzyme
        self.enzymeName = ''

        # Params: Files
        self.fileBg = []
        self.fileExp = []
        self.pathData = 'data'
        self.pathSeqs = os.path.join(self.pathData, 'sequences')
        self.pathFigs = os.path.join(self.pathData, 'figures')
        self.path = ''

        # Params:
        self.seq5Prime = False
        self.seq3Prime = False
        self.subLength = False

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
        self.figEMSquares = figEMSquares
        if figEMSquares:
            self.figSizeEM = (5, 8)  # (width, height)
        else:
            self.figSizeEM = (9.5, 8)
        self.figSize = (9.5, 8)
        self.figSizeMini = (self.figSize[0], 6)
        self.residueLabelType = 2  # 0 = full AA name, 1 = 3-letter code, 2 = 1 letter
        self.labelSizeTitle = 18
        self.labelSizeAxis = 16
        self.labelSizeTicks = 13
        self.lineThickness = 1.5
        self.tickLength = 4
        self.residues = defaultResidues
        self.letters = [residue[2] for residue in self.residues]
        self.colorsAA = NGS.residueColors()

        # Params:
        self. = False
        self. = False
        self. = False

        # Params:
        self. = False
        self. = False
        self. = False

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


    def pressButton(self, message):
        print(f'Received data: {message}')

        return {'key': 'Returned data'}



    def evalDNA(self, data):
        print(f'Received Data:')
        for key, value in data.items():
            print(f'     {key}: {value}')
            if 'file' in key:
                if value:
                    self.loadDNA(value)
        print('\n')

        return {'seq': 'GTGGAACATACCGTGGCGCTGAAACAGAACCGC'}



    def loadDNA(self):
        print()


    def filterAA(self):
        print('Processing Substrates')

        return {'AA': 'VEHTVALKQNR'}



    def filterMotif(self):
        print('Filtering Substrates')

        return {'Motif': 'TVALK'}


