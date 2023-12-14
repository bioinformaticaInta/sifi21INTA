#!/usr/bin/python

import subprocess
from Bio import SeqIO
import os
import json
from collections import Counter
import plotly.graph_objects as go
import plotly.express as px
import plotly.figure_factory as pff
import pandas as pd
from datetime import datetime

import sifi21INTA.generalFunctions as generalFunctions
import sifi21INTA.Sirna as Sirna

class SifiPipeline:
    def __init__(self, bowtieDB, queryFile, outputDir, mode=0, maxGapSize=None, targetsInRegions=False, sirnaSize=21, mismatches=0, accessibilityCheck=True, accessibilityWindow=8, strandCheck=True, endCheck=True, endStabilityTreshold = 1.0, targetSiteAccessibilityTreshold=0.1, terminalCheck=True):

        # Parameters
        self.mode = mode                                                        # Mode, either RNAi design (0) or off-target prediction (1)
        self.sirnaSize = sirnaSize                                              # siRNA size
        self.mismatches = mismatches                                            # Allowed mismatches
        self.strandCheck = strandCheck                                          # Strand selection is enabled or disabled
        self.endCheck = endCheck                                                # End stability selection is enabled or disabled
        self.accessibilityCheck = accessibilityCheck                            # Target site accessibility is enabled or disabled
        self.accessibilityWindow = accessibilityWindow                          # Accessibility window
        self.endStabilityTreshold = endStabilityTreshold                        # End stability treshold
        self.tsAccessibilityTreshold = targetSiteAccessibilityTreshold          # Target site accessibility threshold
        self.terminalCheck = terminalCheck

        self.rnaplfoldLocation = "/usr/local/bin"                               # Rnaplfold path
        self.bowtieLocation = "/usr/bin"                                        # Bowtie path
        self.bowtieDB = bowtieDB                                                # Bowtie DB complete path
        self.allTargets = {}
        self.mainTargets = []
        self.genomeLenghts = {}                                                 # Dictionary with genome information (name:len)
        
        self.outputDir = outputDir                                              # Outpur directory path
        self.queryFile = queryFile                                              # Query sequence in single fasta format complete path
        query = tuple(SeqIO.parse(queryFile, "fasta"))[0]
        self.queryName = str(query.id).replace(" ","_")
        self.querySequence = str(query.seq)
        
        self.targetsInRegions = targetsInRegions                                # True to take targets separated by regions
        self.maxGapSize = maxGapSize                                            # Maximum number of gap to join 2 regions
        if self.targetsInRegions and not self.maxGapSize:
            self.maxGapSize = len(self.querySequence)
        
        self.sirnaFastaFile = self.outputDir+"/"+self.queryName+"."+str(self.sirnaSize)+"sirnas.fasta"
        self.sirnaData = {}
        self.jsonFileName = self.outputDir+"/"+self.queryName+".json"
        self.jsonList = []
        
        # Constants
        self.winsize = 80                                                       # Average the pair probabil. over windows of given size
        self.span = 40                                                          # Set the maximum allowed separation of a base pair to span
        self.temperature = 22                                                   # Temperature for calculation free energy
        self.startPosition = 0
        self.overhang = 2                                                       # siRNA overhang
        self.endNucleotides = 3                                                 # siRNA end nucleotides

    def refGenomeLenghts(self, refGenome):
        if refGenome:
            os.chdir(self.bowtieLocation)
            process = subprocess.Popen(["bowtie-inspect", "-s", self.bowtieDB], stdout=subprocess.PIPE)
            out,err = process.communicate()
            if not process.returncode:
                out = str(out)
                for line in out.split('\\n'):
                    if "Sequence" in line:
                        line = line.split('\\t')
                        self.genomeLenghts[line[1].split(" ")[0]] = int(line[2])
            else:
                raise Exception("Error obtaining genome information")

    def runPipeline(self):
        # Create all siRNAs of size "sirna_size" and save them into multifasta and tab files
        self.createSirnas()
        if self.mode == 0:
            self.designPipeline()
        else:
            self.offTargetPipeline()

    def createSirnas(self):
        #Create all siRNA's of size "self.sirnaSize" of a sequence.
        # Slice over sequence and split into xmers.
        sirnaFastaFile = open(self.sirnaFastaFile, 'w')
        start = self.startPosition
        for start in range(start , len(self.querySequence)-self.sirnaSize+1):
            sequence = self.querySequence[start:start+self.sirnaSize]
            sequence.upper()
            #Add Sirna object to self.sirnaData dict
            self.sirnaData[start+1] = Sirna.Sirna(sequence)
            #Write sirna sequence to multifasta file for bowtie
            sirnaFastaFile.write('>' + str(start+1) + '\n')
            sirnaFastaFile.write(sequence + '\n')
        sirnaFastaFile.close()

    def designPipeline(self):
        # Run BOWTIE against DB
        self.runBowtie()
        
        # Select main targets
        self.selectMainTargets()
        
        # Run RNAplfold
        self.runRnaplfold()
        self.calculateAllSirnasEfficiency()
        
        # Load results to JSON format
        self.createJsonFile()
        
        print("Main target name:",self.mainTargets)
        
        resultTuple = self.getResultSummary()
        self.printResultTable(resultTuple)

        self.efficiencyFigure()

    def offTargetPipeline(self):
        # Run BOWTIE against DB
        self.runBowtie()

        #Plot alignments figure
        self.bowtieFigure()

        # Load results to JSON format
        self.createJsonFile()

        resultTuple = self.getResultSummary()
        self.printResultTable(resultTuple)

    def runBowtie(self):
        """Run BOWTIE alignment."""
        bowtieFile = self.outputDir+"/"+self.queryName+".bowtiehit"
        os.chdir(self.bowtieLocation)
        process = subprocess.Popen(["bowtie", "-a", "-v", str(self.mismatches),  "-y", "-x", self.bowtieDB, "-f", self.sirnaFastaFile, bowtieFile])
        process.wait()
        if os.path.exists(bowtieFile):
            self.loadBowtieData(bowtieFile)

    def loadBowtieData(self, bowtieFileName):
        #Load list with all alignments information contains: [sirnaName, hitName, hitPos, hitStrand, hitMissmatches]
        bowtieAlignments = generalFunctions.bowtieToList(bowtieFileName)
        #Take the targets per regions defined with alignment overlaps
        if self.targetsInRegions:
            #Modify hitName adding the positions of each region
            generalFunctions.addChromosomeRegions(bowtieAlignments, self.maxGapSize, self.sirnaSize)
        #Load bowtie alignments for each Sirna
        for alignment in bowtieAlignments:
            sirnaName = alignment[0]
            hitName = alignment[1]
            self.sirnaData[sirnaName].addBowtieAlignment(*alignment[1:])
            #Count number of hits for each target
            if hitName not in self.allTargets:
                self.allTargets[hitName] = 1
            else:
                self.allTargets[hitName] += 1

    def runRnaplfold(self):
        os.chdir(self.rnaplfoldLocation)
        seq = open(self.queryFile, 'r').read()
        cwd = self.outputDir+"/"+self.queryName+"_RNAplfold"
        os.mkdir(cwd)
        prc_stdout = subprocess.PIPE
        prc = subprocess.Popen(['RNAplfold', '-W', '%d'%self.winsize,'-L', '%d'% self.span, '-u', '%d'%self.sirnaSize, '-T', '%.2f'%self.temperature], stdin=subprocess.PIPE, stdout=prc_stdout, cwd=cwd)
        prc.stdin.write(seq.encode('utf-8'))
        prc.stdin.write('\n'.encode('utf-8'))
        prc.communicate()

        lunpFileName = cwd + '/' + self.queryName + '_lunp'
        if os.path.exists(lunpFileName):
            self.loadRnaplfoldData(lunpFileName)
        
    def loadRnaplfoldData(self, lunpFileName):
        lunpFile = open(lunpFileName)
        for lunpLine in lunpFile:
            if "#" not in lunpLine:
                lunpLine = lunpLine[:-1].split("\t")
                end = int(lunpLine[0])
                if end >= self.sirnaSize:
                    sirnaName = end-self.sirnaSize+1
                    self.sirnaData[sirnaName].addRnaplfoldData(lunpLine[1:])

    def calculateAllSirnasEfficiency(self):
        for sirnaName in self.sirnaData:
            sirna = self.sirnaData[sirnaName]
            if sirnaName < 3:
                sequenceN2 = None
            else:
                sequenceN2 = self.sirnaData[sirnaName-2].getSequence()
            sirna.calculateEfficiency(sequenceN2, self.accessibilityWindow, self.tsAccessibilityTreshold, self.endStabilityTreshold, self.startPosition, self.endNucleotides, self.overhang, self.terminalCheck, self.strandCheck, self.endCheck, self.accessibilityCheck)

    def selectMainTargets(self):
        if self.allTargets:
            allTargetsNumbers = self.bowtieFigure()
            if self.mode == 0:
                mainTargetNumbers = input("Select main targets (comma separated): ").split(",")
                for mainTargetNumber in mainTargetNumbers:
                    if mainTargetNumber in allTargetsNumbers:
                        self.mainTargets.append(allTargetsNumbers[mainTargetNumber])
                    else:
                        raise Exception("Incorrect target number: %s" % (mainTargetNumber))
        else:
            print("No Bowtie alignmets")
    
    def bowtieFigure(self):
        #Obtain all target regions, because we still don't have main targets selected
        #For bowtie figure efficiency not calculed yet, then effCounts contains all zeros
        allTargetsNumbers = {}
        if self.allTargets:
            posX,effCount,allRegions = self.informationForFigure()
            allTargetsNumbers,fig = self.paintTargetRegions(allRegions, "bowtie")
        
            fig.update_layout(
                title='Bowtie alignment regions per target position', # Title
                xaxis_title='Query position', # x-axis name
                yaxis_title='Database target', # y-axis name
                xaxis_tickangle=45,  # Set the x-axis label angle
                showlegend=True,     # Display the legend
            )
            fig.write_html(self.outputDir+"/"+self.queryName+"_mainTargets_selection_plot.html")

        return allTargetsNumbers

    def informationForFigure(self):
        #Return three variables:
        #   posX -> List with positions in a target sequence (X axis)
        #   effCounts -> List with the amount of efficients sirnas at each position (Y axis)
        #   offRegions -> Dictionary with regions with bowtie matches for each target in reference, excluding self.mainTargets

        queryLen = len(self.querySequence)
        posX = [int(i) for i in range(1,queryLen+1)]
        effCount = [0]*queryLen
        offRegions = {}
        
        for sirnaName in sorted(self.sirnaData.keys()):
            startPos = sirnaName-1
            endPos = startPos + (self.sirnaSize-1)
            sirna = self.sirnaData[sirnaName]
            if sirna.getEfficiency():
                for sirnaPos in range(startPos, endPos):
                    effCount[sirnaPos]+=1
            #Obtain offTargets:
            #   Transcriptome reference -> Transcripts
            #   Genome reference -> Chromose regions
            offTargets = sirna.getOffTargets(self.mainTargets)
            for offTarget in offTargets:
                start = startPos
                count = 1
                if offTarget in offRegions:
                    if startPos <= offRegions[offTarget][-1][1]:
                        lastRegion = offRegions[offTarget].pop()
                        start = lastRegion[0]-1
                        count = lastRegion[2]+1
                else:
                    offRegions[offTarget] = []
                offRegions[offTarget].append((start+1 , endPos+1, count))
        return posX,effCount,offRegions 

    def paintTargetRegions(self, regions, mode="bowtie"):
        regionsData = []
        deltaRegionsData = {}
        targetNumbers = {}
        targetNumber = 1
        for target in regions:
            targetNumbers[str(targetNumber)] = target 
            targetName = target
            targetRegion = ""
            if self.targetsInRegions:
                targetName = "_".join(target.split("_")[0:-1])
                targetRegion = target.split("_")[-1]
            deltaRegionsData[target] = []
            for xstart,xend,count in regions[target]:
                regionsData.append({"Target region name":target, 
                                    "Query start":xstart, "Query end":xend, 
                                    "Target ID":targetName, 
                                    "Target region":targetRegion,
                                    "Sirnas in block":count,
                                    "Total sirnas":self.allTargets[target],
                                    "Target number":targetNumber
                                   })
                deltaRegionsData[target].append(xend-xstart)    
            targetNumber += 1
        df = pd.DataFrame(regionsData) 
        hoverData = {"Target region name":False, "Target ID":True, "Target number":False, "Target region":True, "Sirnas in block":False, "Total sirnas":False,"Query start":False,"Query end":False}
        hoverName = None
        if mode=="bowtie":
            hoverData={"Target region name":False, "Target ID":False, "Target number":True, "Target region":True if self.targetsInRegions else False, "Sirnas in block":True, "Total sirnas":True}
            hoverName = "Target ID"

        fig = px.timeline(df, x_start="Query start", x_end="Query end", y="Target region name", color="Target ID", hover_name=hoverName, hover_data=hoverData)#,color_discrete_map=colorByTarget)
        fig.layout.xaxis.type = 'linear'
        for figData in fig.data:
            xData = []
            yNamePre = None
            for yName in figData.y:
                if yName != yNamePre:
                    yNamePre = yName
                    pos = 0
                xData.append(deltaRegionsData[yName][pos])
                pos += 1
            figData.x = xData
        if mode == "bowtie":
            fig.update_layout(yaxis_showticklabels=False)
        return targetNumbers,fig

    def createJsonFile(self):
        jsonFile = open(self.jsonFileName, "w")
        self.dataToJson()
        json.dump(self.jsonList, jsonFile, indent=4)
        jsonFile.close()
    
    def dataToJson(self):
        #Extracts the data from bowtie results and efficiency and put everything into json format.
        for sirnaName in self.sirnaData:
            sirna = self.sirnaData[sirnaName]
            if sirna.bowtieDataTuple():
                for alignment in sirna.bowtieDataTuple():
                    hitName = alignment[0]
                    offTarget = hitName not in self.mainTargets
                    sirnaDict = {"query_name": self.queryName, "sirna_position": sirnaName}
                    sirnaDict.update(sirna.toDict())
                    sirnaDict.update({"is_off_target": offTarget,
                                      "hit_name": hitName, 
                                      "reference_strand_pos": alignment[1],
                                      "strand": alignment[2], 
                                      "mismatches": alignment[3]})
                    self.jsonList.append(sirnaDict)
            # If no alignments, put energy information only in design mode
            elif self.mode == 0:
                sirnaDict = {"query_name": self.queryName, "sirna_position": sirnaName}
                sirnaDict.update(sirna.toDict())
                sirnaDict.update({"is_off_target": False})
                self.jsonList.append(sirnaDict)

    def getResultSummary(self):
        total = 0
        totalWithoutOff = 0
        efficient = 0
        efficientWithoutOff = 0
        for sirnaName in self.sirnaData:
            sirna = self.sirnaData[sirnaName]
            offTargets = sirna.getOffTargets(self.mainTargets)
            effCheck = sirna.getEfficiency()
            if not offTargets:
                totalWithoutOff += 1
            if effCheck:
                efficient += 1
            if effCheck and not offTargets:
                efficientWithoutOff += 1
            total += 1
        if self.mode == 1:
            efficient = efficientWithoutOff = None
        return (total, totalWithoutOff, efficient, efficientWithoutOff)

    def printResultTable(self, resultTuple):
        print("| Total | Total without off targets | Efficients | Efficients without off targets |")
        print("|---------------------------------------------------------------------------------|")
        print("|  %s   |             %s            |    %s      |               %s               |" % (resultTuple))
        print("|---------------------------------------------------------------------------------|")

    def efficiencyFigure(self):
        posX,effCount,offRegions = self.informationForFigure()
        #Add (0,0) initial point to correct figure error when the initial
        #y value is distinct on zero
        posX = [0]+posX
        effCount = [0]+effCount
        if offRegions:
            targetNumbers, fig = self.paintTargetRegions(offRegions, "efficiency")
            fig.add_scatter(x=posX, y=effCount, line_shape='hv', name="Efficient siRNAs", hovertemplate='Query position: %{x}<br>Number of efficient siRNAs: %{y}')
        else:
            df = pd.DataFrame({"Query position":posX , "Number of efficient siRNAs":effCount})
            fig = px.line(df, x="Query position", y="Number of efficient siRNAs", line_shape='hv')#, name="Efficient siRNAs")
        
        fig.update_layout(
            title='siRNAs efficients per target position', # Title
            xaxis_title='Query position', # x-axis name
            yaxis_title='Number of efficient siRNAs', # y-axis name
            xaxis_tickangle=45,  # Set the x-axis label angle
            hovermode="x unified",
            showlegend = True,
            legend_title_text=''
        )
        fig.write_html(self.outputDir+"/"+self.queryName+"_efficiency_plot.html")
