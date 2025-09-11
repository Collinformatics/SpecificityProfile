import base64
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
import gzip
import io
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import os.path
import pandas as pd
import pickle as pk
import queue
import seaborn as sns
import sys
import threading
import time
import warnings



# Generate figures entirely in memory without opening a window
matplotlib.use('Agg')  # Use a non-interactive backend for servers



class WebApp:
    def __init__(self):
        # Params: Dataset
        self.jobParams = {}
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
        self.figures = {}
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
        self.datasetTag = 'Unfiltered'
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
            else:
                import shutil
                # Remove everything inside the directory
                for filename in os.listdir(self.pathFigs):
                    file_path = os.path.join(self.pathFigs, filename)
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)  # delete file or link
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)  # delete subdirectory
                time.sleep(5)
                os.makedirs(self.pathFigs, exist_ok=True)



    @staticmethod
    def createCustomColorMap(colorType):
        colorType = colorType.lower()
        if colorType == 'counts':
            useGreen = True
            if useGreen:
                # Green
                colors = ['#FFFFFF', '#ABFF9B', '#39FF14', '#2E9418', '#2E9418',
                          '#005000']
            else:
                # Orange
                colors = ['white', 'white', '#FF76FA', '#FF50F9', '#FF00F2',
                          '#CA00DF', '#BD16FF']
        elif colorType == 'stdev':
            colors = ['white', 'white', '#FF76FA', '#FF50F9', '#FF00F2', '#CA00DF',
                      '#BD16FF']
        elif colorType == 'word cloud':
            # ,'#F2A900','#2E8B57','black'
            colors = ['#CC5500', '#CC5500', '#F79620', '#FAA338',
                      '#00C01E', '#1D680D', '#003000', 'black']
            # colors = ['#008631','#39E75F','#CC5500','#F79620','black']
        elif colorType == 'em':
            colors = ['navy', 'royalblue', 'dodgerblue', 'lightskyblue', 'white', 'white',
                      'lightcoral', 'red', 'firebrick', 'darkred']
        else:
            print(f'ERROR: Cannot create colormap. '
                  f'Unrecognized colorType parameter: {colorType}\n')
            sys.exit(1)

        # Create colormap
        if len(colors) == 1:
            colorList = [(0, colors[0]), (1, colors[0])]
        else:
            colorList = [(i / (len(colors) - 1), color) for i, color in enumerate(colors)]
        return LinearSegmentedColormap.from_list('custom_colormap', colorList)



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



    @staticmethod
    def encodeFig(fig):
        # Save to a memory buffer instead of disk
        buffer = io.BytesIO()
        plt.savefig(buffer, format='png', bbox_inches='tight')
        buffer.seek(0)

        # Encode as base64 for embedding in HTML
        figBase64 = base64.b64encode(buffer.getvalue()).decode('utf-8')

        return figBase64



    def pressButton(self, message):
        print(f'Received data: {message}')

        return {'key': 'Returned data'}


    def filterAA(self, form):
        self.log()  # Clear the log

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

        # Record job params
        self.jobParams['Enzyme Name'] = self.enzymeName
        self.jobParams['Substrate Length'] = self.seqLength
        self.jobParams['Dataset Tag'] = self.datasetTag




    def initDataStructures(self):
        # Initialize data structures
        self.subsExp = {}
        self.subsBg = {}
        self.xAxisLabel = [f'R{index}' for index in range(1, self.seqLength + 1)]
        self.countsExp = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)
        self.countsBg = pd.DataFrame(0, index=self.AA, columns=self.xAxisLabel)



    def getFilter(self, data):
        if 'filterPos' in data.keys():
            self.filterPos = data['filterPos']
            self.fixAA = {}
            self.exclAA = {}

            # Get filter params
            print('Filter:')
            for key, value in data.items():
                if 'fix' in key:
                    self.fixAA[key] = value
                    print(f'  {key}: {value}')
                if 'excl' in key:
                    self.exclAA[key] = value
                    print(f'  {key}: {value}')

            # Sort filter params
            self.fixAA = dict(sorted(self.fixAA.items()))
            self.exclAA = dict(sorted(self.exclAA.items()))

        self.getDatasetTag()



    def log(self, txt=None):
        if txt is None:
            with open(self.pathLog, 'w'):
                pass
        else:
            with open(self.pathLog, 'a') as log:
                log.write(f'{txt}\n')



    def logInQueue(self, logQueue):
        with open(self.pathLog, 'a') as log:
            while not logQueue.empty():
                log.write(logQueue.get() + '\n')



    def processSubs(self, substrates, datasetType, filteredAA):
        self.log('================================== Substrates '
                 '===================================')
        self.log(f'Substrates: {datasetType}\n')

        # Inspect sequences
        if not filteredAA:
            filteredSubs = {}
            for substrate, count in substrates.items():
                for AA in substrate:
                    if AA not in self.AA:
                        filteredSubs[substrate] = count
            if filteredSubs:
                self.log(f'Filtering Substrates:\n'
                         f'     If a substrate contains an '
                         f'unaccented AA it will be removed.\n'
                         f'     Accepted: {self.AA}\n\n'
                         f'     Removed Substrates:')
                for substrate, count in filteredSubs.items():
                    substrates.pop(substrate, count)
                    self.log(f'          {substrate}: {count}')
                self.log('')

        # Sort data
        substrates = dict(sorted(substrates.items(),
                                 key=lambda item: item[1], reverse=True))

        # Count AAs
        countMatrix = None
        self.log('Substrate Totals:')
        if datasetType == self.datasetTypes['Exp']:
            self.subsExp = substrates
            countMatrix = self.countsExp
            self.countExpTotal = sum(substrates.values())
            self.countExpUnique = len(substrates.keys())
            self.log(f'     Total Substrates: {self.countExpTotal:,}\n'
                     f'    Unique Substrates: {self.countExpUnique:,}\n')

            # Record job params
            self.jobParams['Total Experimental Substrates'] = f'{self.countExpTotal:,}'
            self.jobParams['Unique Experimental Substrates'] = f'{self.countExpUnique:,}'
        elif datasetType == self.datasetTypes['Bg']:
            self.subsBg = substrates
            countMatrix = self.countsBg
            self.countBgTotal = sum(substrates.values())
            self.countBgUnique = len(substrates.keys())
            self.log(f'     Total Substrates: {self.countBgTotal:,}\n'
                     f'    Unique Substrates: {self.countBgUnique:,}\n')

            # Record job params
            self.jobParams['Total Background Substrates'] = f'{self.countBgTotal:,}'
            self.jobParams['Unique Background Substrates'] = f'{self.countBgUnique:,}'
        else:
            self.logErrorFn(function='sampleSize()',
                            msg=f'Unknown dataset type: {datasetType}')

        self.log(f'Top {self.printN:,} {datasetType} Sequences')
        for index, (sub, count) in enumerate(substrates.items()):
            if index >= self.printN:
                break
            self.log(f'     {sub}: {count}')
        self.log('')

        # Save data
        self.saveSubstrates(substrates=substrates,
                            datasetType=datasetType,
                            filteredAA=filteredAA)

        # Count AAs
        self.countAA(substrates=substrates, countMatrix=countMatrix,
                     datasetType=datasetType)


    def logErrorFn(self, function, msg, getStr=False):
        if getStr:
            return (f'\n========================================='
                    f'========================================\n'
                     f'========================================='
                     f'========================================\n\n'
                     f'ERROR: {function}\n'
                     f'{msg}\n\n'
                     f'========================================='
                     f'========================================\n'
                     f'========================================='
                     f'========================================\n')
        else:
            self.log(f'\n========================================='
                     f'========================================\n'
                     f'========================================='
                     f'========================================\n\n'
                     f'ERROR: {function}\n'
                     f'{msg}\n\n'
                     f'========================================='
                     f'========================================\n'
                     f'========================================='
                     f'========================================\n')
            sys.exit(1)



    def evalDNA(self, data):
        self.log()  # Clear the log

        # Get the files
        for key, value in data.items():
            if 'fileExp' in key:
                self.fileExp.append(value)
            elif 'fileBg' in key:
                self.fileBg.append(value)
            else:
                print(key, value)
        print()

        # # Placeholder for files
        self.fileExp = ['data/variantsExp.fastq', 'data/variantsExp2.fastq']
        self.fileBg = ['data/variantsBg.fasta', 'data/variantsBg2.fasta']
        print(f'File Exp: {type(self.fileExp)}\n'
              f'{self.fileExp}\n')
        print(f'File Bg: {type(self.fileBg)}\n'
              f'{self.fileBg}\n\n')


        # Process job parameters
        self.enzymeName = data['enzymeName']
        self.seq5Prime = data['seq5Prime']
        self.seq3Prime = data['seq3Prime']
        self.seqLength = int(data['seqLength'])
        self.minPhred = int(data['minPhred']) if data['minPhred'] != '' else 0


        # Log job params
        self.log(f'================================== Process DNA '
                 f'==================================')
        self.log(f'5\' Sequence: {self.seq5Prime}\n'
                 f'3\' Sequence: {self.seq3Prime}\n'
                 f'Sequence Length: {self.seqLength}\n'
                 f'Min Phred Score: {self.minPhred}')


        # Evaluate job params
        self.initDataStructures()
        self.getFilter(data)
        self.log(f'Dataset Filter: {self.datasetTag}\n\n')


        # Load the data
        threads = []
        queuesExp = []
        queuesExpLog = []
        queuesBg = []
        queuesBgLog = []
        if self.fileExp:
            for file in self.fileExp:
                queueExp = queue.Queue()
                queueLog = queue.Queue()
                queuesExp.append(queueExp)
                queuesExpLog.append(queueLog)
                thread = threading.Thread(
                    target=self.loadDNA,
                    args=(file, self.datasetTypes['Exp'], queueExp, queueLog, True,))
                thread.start()
                threads.append(thread)
        if self.fileBg:
            for file in self.fileBg:
                queueBg = queue.Queue()
                queueLog = queue.Queue()
                queuesBg.append(queueBg)
                queuesBgLog.append(queueLog)
                thread = threading.Thread(
                    target=self.loadDNA,
                    args=(file, self.datasetTypes['Bg'], queueBg, queueLog, True,))
                thread.start()
                threads.append(thread)


        # Wait for all threads to finish
        for thread in threads:
            thread.join()

        # Log the output
        if queuesExpLog:
            for log in queuesExpLog:
                self.logInQueue(log)
        if queuesBgLog:
            for log in queuesBgLog:
                self.logInQueue(log)

        # Get results from queue
        if self.fileExp:
            for queueData in queuesExp:
                substrates = queueData.get()
                for substrate, count in substrates.items():
                    if substrate in self.subsExp.keys():
                        self.subsExp[substrate] += count
                    else:
                        self.subsExp[substrate] = count
        if self.fileBg:
            for queueData in queuesBg:
                substrates = queueData.get()
                for substrate, count in substrates.items():
                    if substrate in self.subsBg.keys():
                        self.subsBg[substrate] += count
                    else:
                        self.subsBg[substrate] = count

        # Make figures
        self.figures = {'exp_counts': False, 'bg_counts': False}
        if self.subsExp:
            # Sort substrates and count AA
            self.processSubs(substrates=self.subsExp,
                             datasetType=self.datasetTypes['Exp'],
                             filteredAA=False)

            # Plot counts
            self.figures['exp_counts'] = (
                self.plotCounts(countedData=self.countsExp,
                                totalCounts=self.countExpTotal,
                                datasetType=self.datasetTypes['Exp']))
        if self.subsBg:
            # Sort substrates and count AA
            self.processSubs(substrates=self.subsBg,
                             datasetType=self.datasetTypes['Bg'],
                             filteredAA=False)

            # Plot counts
            self.figures['bg_counts'] = (
                self.plotCounts(countedData=self.countsBg,
                                totalCounts=self.countBgTotal,
                                datasetType=self.datasetTypes['Bg']))



    def plotCounts(self, countedData, totalCounts, datasetType):
        print(f'Script: {os.path.basename(__file__)}\n')
        print(f"Script: {os.path.abspath(__file__)}\n")


        # Remove commas from string values and convert to float
        countedData = countedData.applymap(lambda x:
                                           float(x.replace(',', ''))
                                           if isinstance(x, str) else x)
        countedData.index = self.AA

        # Create color map
        cMapCustom = self.createCustomColorMap(colorType='Counts')

        # Set figure title
        title = f'\n\n{self.enzymeName}\n{datasetType}\nN={totalCounts:,}'


        # Plot the heatmap with numbers centered inside the squares
        fig, ax = plt.subplots(figsize=self.figSize)
        heatmap = sns.heatmap(countedData, annot=True, fmt=',d', cmap=cMapCustom,
                              cbar=True, linewidths=self.lineThickness-1,
                              linecolor='black', square=False, center=None,
                              annot_kws={'fontweight': 'bold'})
        ax.set_xlabel('Substrate Position', fontsize=self.labelSizeAxis)
        ax.set_ylabel('Residue', fontsize=self.labelSizeAxis)
        ax.set_title(title, fontsize=self.labelSizeTitle, fontweight='bold')
        figBorders = [0.852, 0.075, 0.117, 1]
        plt.subplots_adjust(top=figBorders[0], bottom=figBorders[1],
                            left=figBorders[2], right=figBorders[3])

        # Set the thickness of the figure border
        for _, spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_linewidth(self.lineThickness)

        # Set tick parameters
        ax.tick_params(axis='both', which='major', length=self.tickLength,
                       labelsize=self.labelSizeTicks, width=self.lineThickness)
        ax.tick_params(axis='y', labelrotation=0)

        # Set x-ticks
        xTicks = np.arange(len(countedData.columns)) + 0.5
        ax.set_xticks(xTicks)
        ax.set_xticklabels(countedData.columns)

        # Set y-ticks
        yTicks = np.arange(len(countedData.index)) + 0.5
        ax.set_yticks(yTicks)
        ax.set_yticklabels(countedData.index)


        for _, spine in ax.spines.items():
            spine.set_visible(True)

        # Modify the colorbar
        cbar = heatmap.collections[0].colorbar
        cbar.ax.tick_params(axis='y', which='major', labelsize=self.labelSizeTicks,
                            length=self.tickLength, width=self.lineThickness)
        cbar.outline.set_linewidth(self.lineThickness)
        cbar.outline.set_edgecolor('black')

        # File path
        figName = f'counts - {self.enzymeName} - {datasetType}.png'
        path = os.path.join(self.pathFigs, figName)
        print(f'Saving Fig: {datasetType}\n     {path}\n')

        # Encode the figure
        figBase64 = self.encodeFig(fig)
        with open(path, "wb") as file:
            file.write(base64.b64decode(figBase64))

        # Close the figure to free memory
        plt.close(fig) 
        
        return figName



    def loadDNA(self, path, datasetType, queueData, queueLog, forwardRead):
        translate = True

        # Open the file
        openFn = gzip.open if path.endswith('.gz') else open  # Define open function
        with openFn(path, 'rt') as file:  # 'rt' = read text mode
            if '.fastq' in path or '.fq' in path:
                data = SeqIO.parse(file, 'fastq')
                warnings.simplefilter('ignore', BiopythonWarning)
            elif '.fasta' in path or '.fa' in path:
                data = SeqIO.parse(file, 'fasta')
                warnings.simplefilter('ignore', BiopythonWarning)
            else:
                queueLog.put(self.logErrorFn(
                    function='loadDNA()',
                    msg=f'Unrecognized file\n     {path}',
                    getStr=True))
                translate = False

            # Translate the dna
            if translate:
                substrates = self.translate(data, path, datasetType,
                                            queueLog, forwardRead)
                queueData.put(substrates) # Put the substrates in the queue



    def translate(self, data, fileName, datasetType, queueLog, forwardRead):
        queueLog.put('================================= Translate DNA '
                 '=================================')
        data = list(data)
        substrates = {}
        totalSubsExtracted = 0
        totalSeqsDNA = 0
        self.printN = 10

        queueLog.put(f'File Name: {fileName}')
        if forwardRead:
            queueLog.put('Read Type: Forward Read')
        else:
            queueLog.put('Read Type: Reverse Read')
        queueLog.put(f'  Dataset: {datasetType}')

        # Inspect the file
        useQS = False
        for datapoint in data:
            if 'phred_quality' in datapoint.letter_annotations:
                useQS = True
            break
        queueLog.put(f'  Eval QS: {useQS}\n\n')

        # Inspect the datasetType parameter
        if (datasetType != self.datasetTypes['Exp']
                and datasetType != self.datasetTypes['Bg']):
            self.logErrorFn(function='translate()',
                            msg=f'Unknown dataset type: {datasetType}')


        def extractionEfficiency(fullSet=False):
            perExtracted = (totalSubsExtracted / totalSeqsDNA) * 100
            if fullSet:
                queueLog.put('- All Sequences')
            else:
                queueLog.put(f'\nExtraction Efficiency: {datasetType}\n'
                         f'- First {self.printN} Sequences')
            queueLog.put(f'     Evaluated DNA Sequences: {totalSeqsDNA:,}\n'
                         f'        Extracted Substrates: {totalSubsExtracted:,}\n'
                         f'       Extraction Efficiency: {round(perExtracted, 3)} %\n')


        # Translate DNA - Sample Set
        if useQS:
            for index, datapoint in enumerate(data):
                if totalSubsExtracted >= self.printN: # Exit loop
                    break

                # Process datapoint
                totalSeqsDNA += 1
                dna = str(datapoint.seq)
                queueLog.put(f'DNA Seq: {dna}')


                # Inspect full dna seq
                if self.seq5Prime in dna and self.seq3Prime in dna:
                    qs = datapoint.letter_annotations['phred_quality']

                    # Find: Substrate indices
                    start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                    end = dna.find(self.seq3Prime)

                    # Extract substrate dna seq
                    substrateDNA = dna[start:end].strip()
                    queueLog.put(f'Sub DNA: {substrateDNA}')
                    if len(substrateDNA) == self.seqLength * 3:
                        # Express substrate
                        substrate = str(Seq.translate(substrateDNA))
                        queueLog.put(f'Sub Seq: {substrate}')

                        # Inspect substrate seq: PRINT ONLY
                        if 'X' not in substrate and '*' not in substrate:
                            qs = qs[start:end]
                            queueLog.put(f'     QS: {qs}')
                            if all(score >= self.minPhred for score in qs):
                                queueLog.put(f'Keep Substrate\n')
                                if substrate in substrates.keys():
                                    substrates[substrate] += 1
                                else:
                                    substrates[substrate] = 1
                                totalSubsExtracted += 1
                            else:
                                queueLog.put('')
        else:
            for index, datapoint in enumerate(data):

                if totalSubsExtracted >= self.printN:
                    break

                # Process datapoint
                totalSeqsDNA += 1
                dna = str(datapoint.seq)
                queueLog.put(f'DNA Seq: {dna}')

                # Inspect full dna seq
                if self.seq5Prime in dna and self.seq3Prime in dna:

                    # Find: Substrate indices
                    start = dna.find(self.seq5Prime) + len(self.seq5Prime)
                    end = dna.find(self.seq3Prime)

                    # Extract substrate dna seq
                    substrateDNA = dna[start:end].strip()
                    queueLog.put(f'Sub DNA: {substrateDNA}')

                    if len(substrateDNA) == self.seqLength * 3:
                        # Express substrate
                        substrate = str(Seq.translate(substrateDNA))
                        queueLog.put(f'Sub Seq: {substrate}')

                        # Inspect substrate seq: PRINT ONLY
                        if 'X' not in substrate and '*' not in substrate:
                            queueLog.put(f'Keep Substrate\n')
                            if substrate in substrates.keys():
                                substrates[substrate] += 1
                            else:
                                substrates[substrate] = 1
                            totalSubsExtracted += 1
                        else:
                            queueLog.put('')
        extractionEfficiency() # Evaluate data quality


        # Translate DNA - Full Set
        substrates = {}
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
        queueLog.put('')

        return substrates



    def saveSubstrates(self, substrates, datasetType, filteredAA):
        saveTag = None
        if filteredAA:
            if self.datasetTag is None:
                print(f'Dont save, dataset tag: {self.datasetTag}\n')
                sys.exit()

            print(f'Saving Substrates: {datasetType}\n')
            if datasetType == self.datasetTypes['Exp']:
                saveTag = self.saveTagExp['subs']
            elif datasetType == self.datasetTypes['Bg']:
                saveTag = self.saveTagBg['subs']
            else:
                self.logErrorFn(function='saveSubstrates()',
                                msg=f'Unknown dataset type: {datasetType}')
        else:
            if datasetType == self.datasetTypes['Exp']:
                saveTag = self.saveTagExp['subsRaw']
            elif datasetType == self.datasetTypes['Bg']:
                saveTag = self.saveTagBg['subsRaw']
            else:
                self.logErrorFn(function='saveSubstrates()',
                                msg=f'Unknown dataset type: {datasetType}')

        # Save the substrates
        path = os.path.join(self.pathSeqs, saveTag)
        self.log(f'Saving Substrates: {datasetType}\n     {path}\n\n')
        with open(path, 'wb') as file:
            pk.dump(substrates, file)



    def countAA(self, substrates, countMatrix, datasetType):
        self.log('=================================== Count AA '
                 '====================================')
        self.log(f'Dataset: {datasetType}\n'
                 f'Unique Substrates: {len(substrates.keys())}')
        totalCounts = 0
        for substrate, count in substrates.items():
            totalCounts += count
            for index, AA in enumerate(substrate):
                countMatrix.loc[AA, f'R{index + 1}'] += count
        self.log(f'Unique Substrates: {len(substrates.keys())}\n'
                 f'Total Counts: {totalCounts}\n\n{countMatrix}\n\n')


