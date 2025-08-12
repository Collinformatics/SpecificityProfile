from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
import gzip
import os.path
import pandas as pd
import queue
import sys
import threading
import warnings



class WebApp:
    def __init__(self):
        # Params: Dataset
        self.enzymeName = ''
        self.seqLength = False
        self.minCounts = 1
        self.printN = 10
        self.xAxisLabel = False
        self.subsExp = {}
        self.countsExp = 'Initialize me'
        self.countExpTotal = 0
        self.countExpUnique = 0
        self.saveTagExp = {}
        self.subsBg = {}
        self.countsBg = 'Initialize me'
        self.countBgTotal = 0
        self.countBgUnique = 0
        self.saveTagBg = {}
        self.saveTagFig = {}
        self.fixMotif = False
        self.datasetTypes = {'Exp': 'Experimental',
                             'Bg': 'Background'}

        # Params: Files
        self.queueLog = queue.Queue()
        self.fileExp = []
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
        self.figSize = (9.5, 8)
        self.figSizeMini = (self.figSize[0], 6)
        self.residueLabelType = 2  # 0 = full AA name, 1 = 3-letter code, 2 = 1 letter
        self.labelSizeTitle = 18
        self.labelSizeAxis = 16
        self.labelSizeTicks = 13
        self.lineThickness = 1.5
        self.tickLength = 4
        self.colorsAA = self.residueColors()
        self.AA = list(self.colorsAA.keys())

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


    def filterAA(self, form):
        self.log()  # Clear the log

        # Process job parameters
        self.enzymeName = form['enzymeName']
        self.seqLength = int(form['seqLength'])
        print('Processing Substrates')

        # Evaluate job params
        self.getFilter(form)
        
        return {'AA': 'VEHTVALKQNR'}


    def filterMotif(self, form):
        self.log()  # Clear the log

        # Process job parameters
        self.enzymeName = form['enzymeName']
        self.seqLength = int(form['seqLength'])
        print('Filtering Substrates')
        self.fixMotif = True

        # Evaluate job params
        self.getFilter(form)

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


        # Initialize: Save tags
        self.saveTagExp = {
            'subsRaw': f'{self.enzymeName} - Subs Exp - '
                    f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'countsRaw': f'{self.enzymeName} - AA Counts Exp - '
                      f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'subs': f'{self.enzymeName} - Subs Exp - {self.datasetTag} - '
                    f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'counts': f'{self.enzymeName} - AA Counts Exp - {self.datasetTag} - '
                      f'MinCounts {self.minCounts} - {self.seqLength} AA'
        }
        self.saveTagBg = {
            'subsRaw': f'{self.enzymeName} - Subs Exp - '
                   f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'countsRaw': f'{self.enzymeName} - AA Counts Bg - '
                      f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'subs': f'{self.enzymeName} - Subs Bg - {self.datasetTag} - '
                    f'MinCounts {self.minCounts} - {self.seqLength} AA',
            'counts': f'{self.enzymeName} - AA Counts Bg - {self.datasetTag} - '
                      f'MinCounts {self.minCounts} - {self.seqLength} AA',
        }
        if self.fixMotif:
            for tag, path in self.saveTagExp:
                self.saveTagExp[tag] = (path.replace(f'{self.enzymeName}',
                                                     f'{self.enzymeName} - Motif'))


        self.saveTagFig = (f'{self.enzymeName} - Fig - {self.datasetTag} - '
                           f'Min Counts {self.minCounts} - {self.seqLength} AA')

        # Initialize data structures
        self.xAxisLabel = [f'R{index}' for index in range(1, self.seqLength + 1)]
        self.countsExp = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)
        self.countsBg = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)



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
        # else:
        #     with open(self.pathLog, 'a') as log:
        #         log.write(f'{txt}\n')


    def logInQueue(self, logQueue):
        with open(self.pathLog, 'a') as log:
            while not logQueue.empty():
                log.write(logQueue.get() + '\n')



    def logSubs(self, substrates, datasetType):
        self.log('================================== Substrates '
                 '===================================')
        self.log(f'Substrates: {datasetType}')
        if datasetType == self.datasetTypes['Exp']:
            self.countExpTotal = sum(substrates.values())
            self.countExpUnique = len(substrates.keys())
            self.log(f'     Total Substrates: {self.countExpTotal:,}\n'
                     f'    Unique Substrates: {self.countExpUnique:,}\n')
        elif datasetType == self.datasetTypes['Bg']:
            self.countBgTotal = sum(substrates.values())
            self.countBgUnique = len(substrates.keys())
            self.log(f'     Total Substrates: {self.countBgTotal:,}\n'
                     f'    Unique Substrates: {self.countBgUnique:,}\n')
        else:
            self.logError(function='sampleSize()',
                          msg=f'Unknown dataset type: {datasetType}')

        self.log(f'Top {self.printN:,} {datasetType} Sequences')
        for index, (sub, count) in enumerate(substrates.items()):
            if index >= self.printN:
                self.log('')
                break
            self.log(f'     {sub}: {count}')



    def logError(self, function, msg):
        self.log(f'\n========================================='
                 f'========================================\n'
                 f'========================================='
                 f'========================================\n\n'
                 f'ERROR: {function}\n'
                 f'{msg}\n\n'
                 f'========================================='
                 f'========================================\n'
                 f'========================================='
                 f'========================================\n'
                 )
        sys.exit(1)



    def evalDNA(self, form):
        self.log()  # Clear the log

        # Process job parameters
        self.enzymeName = form['enzymeName']
        self.fileExp = form['fileExp']
        self.fileBg = form['fileBg']
        self.seq5Prime = form['seq5Prime']
        self.seq3Prime = form['seq3Prime']
        self.seqLength = int(form['seqLength'])
        self.minPhred = int(form['minPhred']) if form['minPhred'] != '' else 0


        # Log job params
        self.log(f'================================== Process DNA '
                 f'==================================')
        self.log(f'5\' Sequence: {self.seq5Prime}\n'
                 f'3\' Sequence: {self.seq3Prime}\n'
                 f'Sequence Length: {self.seqLength}\n'
                 f'Min Phred Score: {self.minPhred}')


        # Evaluate job params
        self.getFilter(form)
        self.log(f'Dataset Filter: {self.datasetTag}\n\n')

        # Params: Tmp
        # Update when using files

        # if type(file) == 'fastq':
        #     useQS = True

        self.fileExp = ['data/variantsExp.fastq', 'data/variantsExp2.fastq']
        self.fileBg = ['data/variantsBg.fasta', 'data/variantsBg2.fasta']
        self.fileBg = False

        # Load the data
        queueExp = queue.Queue()
        queueBg = queue.Queue()
        threads = []
        logQueues = []
        if self.fileExp:
            # self.loadDNA(path=self.fileExp,
            #              datasetType=self.datasetTypes['Exp'],
            #              forwardRead=True)
            for file in self.fileExp:
                queueLog = queue.Queue()
                logQueues.append((f'exp_{file}.log', queueLog))
                thread = threading.Thread(
                    target=self.loadDNA,
                    args=(file, self.datasetTypes['Exp'], queueExp, queueLog, True,))
                thread.start()
                threads.append(thread)
        if self.fileBg:
            # self.loadDNA(path=self.fileBg,
            #              datasetType=self.datasetTypes['Bg'],
            #              forwardRead=True)
            for file in self.fileBg:
                thread = threading.Thread(
                    target=self.loadDNA,
                    args=(file, self.datasetTypes['Bg'], queueBg, True,))
            thread.start()
            threads.append(thread)

        # Wait for all threads to finish
        for thread in threads:
            thread.join()
        print(f'Here')

        # Get results from queue
        print(queueExp)
        if self.fileExp:
            print('Start A')
            resultsExp = queueExp.get()
            print('Done A')
        if self.fileBg:
            print('Start B')
            resultsBg = queueBg.get()
            print('Done B')
        print('There')

        if logQueues:
            for logName, queueLog in logQueues:
                print(f'Log: {logName}, {type(queueLog)}')
                self.logInQueue(queueLog)

        print(f'\nResults:')
        for x in resultsExp:
            print(x)
        sys.exit()

        print(f'\nResults: {resultsExp}\n\n')
        print(f'\nResults: {resultsBg}\n\n')


        # Sort data
        if self.subsExp:
            self.subsExp = dict(sorted(self.subsExp.items(),
                                       key=lambda item: item[1], reverse=True))
            self.logSubs(substrates=self.subsExp, datasetType=self.datasetTypes['Exp'])
            self.saveSubstrates(substrates=self.subsExp,
                                datasetType=self.datasetTypes['Exp'], rawSubs=True)
        if self.subsBg:
            self.subsBg = dict(sorted(self.subsBg.items(),
                                      key=lambda item: item[1], reverse=True))
            self.logSubs(substrates=self.subsBg, datasetType=self.datasetTypes['Bg'])
            self.saveSubstrates(substrates=self.subsBg,
                                datasetType=self.datasetTypes['Bg'], rawSubs=True)


        # Count AAs
        if self.fileExp:
            self.countAA(substrates=self.subsExp, countMatrix=self.countsExp,
                         datasetType=self.datasetTypes['Exp'], fromDNA=True)
        if self.fileBg:
            self.countAA(substrates=self.subsBg, countMatrix=self.countsBg,
                         datasetType=self.datasetTypes['Bg'], fromDNA=True)

        return {'seq': 'GTGGAACATACCGTGGCGCTGAAACAGAACCGC'}



    def loadDNA(self, path, datasetType, queueData, queueLog, forwardRead):
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
                queueLog.put(f'ERROR: Unrecognized file\n     {path}\n\n')
                # self.log(f'ERROR: Unrecognized file\n     {path}\n\n')
            queueLog.put(f'Translate: {data}')
            # self.log(f'Translate: {data}')

            # Translate the dna
            substrates = self.translate(data, datasetType, queueLog, forwardRead)
            queueData.put(substrates) # Put substrates in the queue



    def translate(self, data, datasetType, queueLog, forwardRead):
        self.log('================================= Translate DNA '
                 '=================================')
        queueLog.put('================================= Translate DNA '
                 '=================================')
        data = list(data)

        substrates = {}
        totalSubsExtracted = 0
        totalSeqsDNA = 0
        self.printN = 10
        if forwardRead:
            self.log('Translating: Forward Read')
            queueLog.put('Translating: Forward Read')
        else:
            self.log('Translating: Reverse Read')
            queueLog.put()
        self.log(f'    Dataset: {datasetType}')
        queueLog.put(f'    Dataset: {datasetType}')

        # Inspect the file
        useQS = False
        for datapoint in data:
            if 'phred_quality' in datapoint.letter_annotations:
                useQS = True
            break
        self.log(f' Inspect QS: {useQS}\n\n')
        queueLog.put(f' Inspect QS: {useQS}\n\n')

        # Inspect the datasetType parameter
        if (datasetType != self.datasetTypes['Exp']
                and datasetType != self.datasetTypes['Bg']):
            self.logError(function='translate()',
                          msg=f'Unknown dataset type: {datasetType}')


        def extractionEfficiency(fullSet=False):
            perExtracted = (totalSubsExtracted / totalSeqsDNA) * 100
            if fullSet:
                self.log('- All Sequences')
                queueLog.put('- All Sequences')
            else:
                self.log(f'\nExtraction Efficiency: {datasetType}\n'
                         f'- First {self.printN} Sequences')
                queueLog.put(f'\nExtraction Efficiency: {datasetType}\n'
                         f'- First {self.printN} Sequences')
            self.log(f'     Evaluated DNA Sequences: {totalSeqsDNA:,}\n'
                     f'        Extracted Substrates: {totalSubsExtracted:,}\n'
                     f'       Extraction Efficiency: {round(perExtracted, 3)} %\n')
            queueLog.put(f'     Evaluated DNA Sequences: {totalSeqsDNA:,}\n'
                         f'        Extracted Substrates: {totalSubsExtracted:,}\n'
                         f'       Extraction Efficiency: {round(perExtracted, 3)} %\n')


        # Translate DNA - Sample Set
        print(f'QS: {useQS}\n')
        if useQS:
            print(f'Start: {type(data)}\n')
            for index, datapoint in enumerate(data):
                if totalSubsExtracted >= self.printN: # Exit loop
                    break

                # Process datapoint
                totalSeqsDNA += 1
                dna = str(datapoint.seq)
                self.log(f'DNA Seq: {dna}')
                print(f'DNA Seq: {dna}')

                # Inspect full dna seq
                if self.seq5Prime in dna and self.seq3Prime in dna:
                    qs = datapoint.letter_annotations['phred_quality']

                    # Find: Substrate indices
                    start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                    end = dna.find(self.seq3Prime)

                    # Extract substrate dna seq
                    substrateDNA = dna[start:end].strip()
                    self.log(f'Sub Seq: {substrateDNA}')
                    if len(substrateDNA) == self.seqLength * 3:
                        # Express substrate
                        substrate = str(Seq.translate(substrateDNA))
                        self.log(f'    Sub: {substrate}')
                        print(f'    Sub: {substrate}')

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
                                queueLog.put('')
        else:
            for index, datapoint in enumerate(data):

                if totalSubsExtracted >= self.printN:
                    break

                # Process datapoint
                totalSeqsDNA += 1
                dna = str(datapoint.seq)
                self.log(f'DNA Seq: {dna}')

                # Inspect full dna seq
                if self.seq5Prime in dna and self.seq3Prime in dna:

                    # Find: Substrate indices
                    start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                    end = dna.find(self.seq3Prime)

                    # Extract substrate dna seq
                    substrateDNA = dna[start:end].strip()
                    self.log(f'Sub Seq: {substrateDNA}')

                    if len(substrateDNA) == self.seqLength * 3:
                        # Express substrate
                        substrate = str(Seq.translate(substrateDNA))
                        self.log(f'    Sub: {substrate}')

                        # Inspect substrate seq: PRINT ONLY
                        if 'X' not in substrate and '*' not in substrate:
                            self.log(f'Keep Substrate\n')
                            if substrate in substrates.keys():
                                substrates[substrate] += 1
                            else:
                                substrates[substrate] = 1
                            totalSubsExtracted += 1
                        else:
                            self.log('')
        extractionEfficiency() # Evaluate data quality


        # Translate DNA - Full Set
        totalSeqsDNA = 0
        totalSubsExtracted = 0
        if useQS:
            for index, datapoint in enumerate(data):
                # Process datapoint
                totalSeqsDNA += 1
                dna = str(datapoint.seq)

                # Inspect full dna seq
                if self.seq5Prime in dna and self.seq3Prime in dna:
                    qs = datapoint.letter_annotations['phred_quality']

                    # Find: Substrate indices
                    start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                    end = dna.find(self.seq3Prime)

                    # Extract substrate dna seq
                    substrateDNA = dna[start:end].strip()
                    if len(substrateDNA) == self.seqLength * 3:
                        # Express substrate
                        substrate = str(Seq.translate(substrateDNA))

                        # Inspect substrate seq: PRINT ONLY
                        if 'X' not in substrate and '*' not in substrate:
                            qs = qs[start:end]
                            if all(score >= self.minPhred for score in qs):
                                if substrate in substrates.keys():
                                    substrates[substrate] += 1
                                else:
                                    substrates[substrate] = 1
                                totalSubsExtracted += 1
        else:
            for index, datapoint in enumerate(data):
                # Process datapoint
                totalSeqsDNA += 1
                dna = str(datapoint.seq)

                # Inspect full dna seq
                if self.seq5Prime in dna and self.seq3Prime in dna:

                    # Find: Substrate indices
                    start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                    end = dna.find(self.seq3Prime)

                    # Extract substrate dna seq
                    substrateDNA = dna[start:end].strip()
                    if len(substrateDNA) == self.seqLength * 3:
                        # Express substrate
                        substrate = str(Seq.translate(substrateDNA))

                        # Inspect substrate seq: PRINT ONLY
                        if 'X' not in substrate and '*' not in substrate:
                            if substrate in substrates.keys():
                                substrates[substrate] += 1
                            else:
                                substrates[substrate] = 1
                            totalSubsExtracted += 1
        extractionEfficiency(fullSet=True)  # Evaluate data quality
        self.log('')
        print('end:', totalSeqsDNA, len(data))
        queueLog.put(f'End:, {totalSeqsDNA}, {len(data)}')

        return substrates



    def saveSubstrates(self, substrates, datasetType, rawSubs=False):
        path = None
        if rawSubs:
            if datasetType == self.datasetTypes['Exp']:
                path = self.saveTagExp['subsRaw']
            elif datasetType == self.datasetTypes['Bg']:
                path = self.saveTagBg['subsRaw']
            else:
                self.logError(function='saveSubstrates()',
                              msg=f'Unknown dataset type: {datasetType}')
        elif self.datasetTag is not None:
            print(f'Saving Substrates: {datasetType}\n')
            if datasetType == self.datasetTypes['Exp']:
                path = self.saveTagExp['subs']
            elif datasetType == self.datasetTypes['Bg']:
                path = self.saveTagBg['subs']
            else:
                self.logError(function='saveSubstrates()',
                              msg=f'Unknown dataset type: {datasetType}')
        else:
            print(f'Dont save, dataset tag: {self.datasetTag}\n')

        self.log(f'Saving Substrates: {datasetType}\n     {path}\n\n')



    def countAA(self, substrates, countMatrix, datasetType, fromDNA=False):
        self.log('=================================== Count AA '
                 '====================================')
        self.log(f'Dataset: {datasetType}\n'
                 f'Unique Substrates: {len(substrates.keys())}')

        if fromDNA:
            for substrate, count in substrates.items():
                countSub = True
                for AA in substrate:
                    if AA not in self.AA:
                        countSub = False
                        self.log(f'Warning: An AA ({AA}) in substrate ({substrate}) '
                                 f'is not an accepted AAs\n'
                                 f'     Accepted: {self.AA}\n')
                        break
                if countSub:
                    for index, AA in enumerate(substrate):
                        countMatrix.loc[AA, f'R{index + 1}'] += count
        self.log(f'Total {countMatrix}\n\n')

