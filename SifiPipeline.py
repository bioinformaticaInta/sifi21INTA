#!/usr/bin/python

import subprocess
from Bio import SeqIO
import os
import shutil
import json
import plotly.express as px
import pandas as pd
import glob

import sifi21INTA.generalFunctions as generalFunctions
import sifi21INTA.Sirna as Sirna
import sifi21INTA.TargetSequence as TargetSequence

class SifiPipeline:
    def __init__(self, bowtieDB, queryFile, outputDir, mode=0, maxGapSize=None, genomicReference=True, gffFile=None, targetsInRegions=False, sirnaSize=21, mismatches=0, accessibilityCheck=True, accessibilityWindow=8, strandCheck=True, endCheck=True, endStabilityTreshold = 1.0, targetSiteAccessibilityTreshold=0.1, terminalCheck=True):

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

        self.bowtieDB = bowtieDB                                                # Bowtie DB path
        self.allTargets = {}
        self.mainTargets = []
        
        self.outputDir = outputDir                                              # Outpur directory path
        self.queryFile = queryFile                                              # Query sequence in single fasta format path
        query = tuple(SeqIO.parse(queryFile, "fasta"))[0]
        self.queryName = str(query.id).replace(" ","_")
        self.querySequence = str(query.seq)
        
        self.gffFile = gffFile                                                  # GFF annotation file path
        self.genomicReference = genomicReference                                # True if the reference is a genome
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

    def runPipeline(self):
        # Create all siRNAs of size "sirna_size" and save them into multifasta and tab files
        self.createSirnas()
        if self.mode == 0:
            self.designPipeline()
        else:
            self.offTargetPipeline()

    def designPipeline(self):
        # Run BOWTIE against DB
        print()
        print("#####################################################")
        print("Step 1: Run BOWTIE alignment")
        print("Selected Database: "+self.bowtieDB)
        print("Query file: "+self.queryFile)
        print("#####################################################")
        print()
        self.runBowtie()

        # Select main targets
        print()
        print("#####################################################")
        print("Step 2: Select main target for siRNAs prediction")
        print("#####################################################")
        print()
        self.selectMainTargets()
        print("Main target selected:",self.mainTargets)
        print()

        # Run RNAplfold
        print()
        print("#####################################################")
        print("Step 3: Calculate efficiency")
        print("Using default method and parameters")
        print("#####################################################")
        print()
        self.runRnaplfold()
        self.calculateAllSirnasEfficiency()
        
        print()
        print("#####################################################")
        print("Step 4: Results")
        print("#####################################################")
        print()
        
        # Load results to JSON format
        self.createJsonFile()
        print("Created JSON output file: " + self.jsonFileName)
        print()
        
        print("Summary table:")
        resultTuple = self.getResultSummary()
        self.printResultTable(resultTuple)
        print()

        figureName = self.efficiencyFigure()
        print("Efficiency figure created: " + figureName)
        print()

    def offTargetPipeline(self):
        # Run BOWTIE against DB
        print()
        print("#####################################################")
        print("Step 1: Run BOWTIE alignment")
        print("Selected Database: "+self.bowtieDB)
        print("Query file: "+self.queryFile)
        print("#####################################################")
        print()
        self.runBowtie()

        #Plot alignments figure
        allTargetsNumbers,figureName = self.bowtieFigure()
        print()
        print("#####################################################")
        print("Step 2: Plot alignments figure")
        print("Alignment figure created: " + figureName)
        print("#####################################################")
        print()

        print()
        print("#####################################################")
        print("Step 3: Results")
        print("#####################################################")
        print()
        
        # Load results to JSON format
        self.createJsonFile()
        print("Created JSON output file: " + self.jsonFileName)
        print()
        
        print("Summary table:")
        resultTuple = self.getResultSummary()
        self.printResultTable(resultTuple)
        print()

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

    def runBowtie(self):
        """Run BOWTIE alignment."""
        bowtieFile = self.outputDir+"/"+self.queryName+".bowtiehit"
        process = subprocess.Popen(["bowtie", "-a", "-v", str(self.mismatches),  "-y", "-x", self.bowtieDB, "-f", self.sirnaFastaFile, bowtieFile])
        process.wait()
        if os.path.exists(bowtieFile):
            self.loadBowtieData(bowtieFile)

    def loadBowtieData(self, bowtieFileName):
        #Load list with all alignments information contains: [sirnaName, hitTarget, hitPos, hitStrand, hitMissmatches]
        bowtieAlignments, targetsRegions = generalFunctions.bowtieToList(bowtieFileName)
        #Take the targets per regions defined with alignment overlaps
        if bowtieAlignments:
            if self.targetsInRegions:
                #Add the positions of the region for each alignment in the reference
                targetsRegions = generalFunctions.addTargetsRegions(bowtieAlignments, self.maxGapSize, self.sirnaSize)
            #Add annotation from gff file
            if self.gffFile:
                generalFunctions.addAnnotations(self.gffFile, bowtieAlignments, self.genomicReference, self.targetsInRegions, targetsRegions)
            #Load bowtie alignments for each Sirna
            countedSirnas = set()
            for alignment in bowtieAlignments:
                sirnaName = alignment[0]
                hitTarget = alignment[1]
                self.sirnaData[sirnaName].addBowtieAlignment(*alignment[1:])
                #Count number of hits for each target only one time for each
                if hitTarget not in self.allTargets:
                    self.allTargets[hitTarget] = 1
                elif (sirnaName,hitTarget) not in countedSirnas:
                    self.allTargets[hitTarget] += 1
                countedSirnas.add((sirnaName,hitTarget))

    def runRnaplfold(self):
        seq = open(self.queryFile, 'r').read()
        cwd = self.outputDir+"/"+self.queryName+"_RNAplfold"
        os.mkdir(cwd)
        prc_stdout = subprocess.PIPE
        prc = subprocess.Popen(['RNAplfold', '-W', '%d'%self.winsize,'-L', '%d'% self.span, '-u', '%d'%self.sirnaSize, '-T', '%.2f'%self.temperature], stdin=subprocess.PIPE, stdout=prc_stdout, cwd=cwd)
        prc.stdin.write(seq.encode('utf-8'))
        prc.stdin.write('\n'.encode('utf-8'))
        prc.communicate()

        lunpFileName = glob.glob(cwd+'/*lunp')
        if lunpFileName:
            self.loadRnaplfoldData(lunpFileName[0])
        shutil.rmtree(cwd)

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
            allTargetsNumbers,figureName = self.bowtieFigure()
            if self.mode == 0:
                mainTargetNumbers = input("Select main targets observing the figure: "+ figureName +" (insert several numbers comma separated): ").split(",")
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
            figureName = self.outputDir+"/"+self.queryName+"_mainTargets_selection_plot.html"
            fig.write_html(figureName)
        return allTargetsNumbers,figureName

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
        figureName = self.outputDir+"/"+self.queryName+"_efficiency_plot.html"
        fig.write_html(figureName)
        return figureName
    
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
            startQueryPos = sirnaName-1
            endQueryPos = startQueryPos + (self.sirnaSize-1)
            sirna = self.sirnaData[sirnaName]
            if sirna.getEfficiency():
                for sirnaPos in range(startQueryPos, endQueryPos):
                    effCount[sirnaPos]+=1
            #Obtain offTargets:
            #   Transcriptome reference -> Transcripts
            #   Genome reference -> Chromosome regions
            offTargets = sirna.getOffTargets(self.mainTargets)
            for offTarget,startInOffTarget,strandInOffTarget,missmtcInOffTarget in offTargets:
                startQuery = startQueryPos
                endQuery = endQueryPos
                endInOffTarget = startInOffTarget + (self.sirnaSize-1)
                count = 1
                #Identify the block where the alignment from in the off Target sequence
                if offTarget in offRegions:
                    pos = 0
                    found = False
                    while pos < len(offRegions[offTarget]) and not found:
                        #The alingment belongs to this block, then extend the extremes
                        #start extreme in negative strand or end extreme in positive strand 
                        offRegion = offRegions[offTarget][pos]
                        if strandInOffTarget == "+":
                            if offRegion[0] <= startQueryPos <= offRegion[1] and offRegion[2] <= startInOffTarget <= offRegion[3] and strandInOffTarget == offRegion[5]:
                                found = True
                                startInOffTarget = offRegion[2]
                        else:
                            if offRegion[0] <= endQueryPos <= offRegion[1] and offRegion[2] <= endInOffTarget <= offRegion[3] and strandInOffTarget == offRegion[5]:
                                found = True
                                endInOffTarget = offRegion[3]
                        if found:
                            offRegion = offRegions[offTarget].pop(pos)
                            startQuery = offRegion[0]-1
                            count = offRegion[4]+1
                        pos += 1
                else:
                    offRegions[offTarget] = []
                offRegions[offTarget].append((startQuery+1 , endQuery+1, startInOffTarget, endInOffTarget, count, strandInOffTarget))
        return posX,effCount,offRegions 

    def paintTargetRegions(self, regions, mode="bowtie"):
        regionsData = []
        deltaRegionsData = {}
        targetNumbers = {}
        targetNumber = 1
        gffData = []
        for target in regions:
            targetNumbers[str(targetNumber)] = target 
            deltaRegionsData[target.getNameWithRegion()] = []
            for qstart,qend,tstart,tend,count,strand in regions[target]:
                regionsData.append({"Target region name":target.getNameWithRegion(), 
                                    "Target number":targetNumber,
                                    "Start position in query":qstart, "End position in query":qend, 
                                    "Target ID":target.getName(), 
                                    "Target region":"-".join(map(str,target.getRegion())),
                                    "Target region length":target.getLength(),
                                    "Alignment strand":strand,
                                    "Start block position in target":tstart, "End block position in target":tend,
                                    "Number of sirnas in block":count,
                                    "Total sirnas in target region" if self.targetsInRegions else "Total sirnas in target":self.allTargets[target],
                                    "Total blocks in target region" if self.targetsInRegions else "Total blocks in target":len(regions[target]),
                                    "Target region annotation" if self.targetsInRegions else "Target annotation":"<br>".join(target.getAnnotation())
                                   })
                deltaRegionsData[target.getNameWithRegion()].append(qend-qstart)    
            targetNumber += 1
        df = pd.DataFrame(regionsData) 
        hoverData = {"Target region name":False, "Target ID":True, "Target number":False, "Target region":True if self.targetsInRegions else False, "Target region length":False, "Alignment strand":False, "Start block position in target":False, "End block position in target":False, "Number of sirnas in block":False, "Total sirnas in target region" if self.targetsInRegions else "Total sirnas in target":False, "Total blocks in target region" if self.targetsInRegions else "Total blocks in target":False, "Start position in query":False,"End position in query":False, "Target region annotation" if self.targetsInRegions else "Target annotation":True if self.gffFile else False}
        hoverName = None
        if mode=="bowtie":
            hoverData = {"Target region name":False, "Target ID":False, "Target number":True, "Target region":True if self.targetsInRegions else False, "Target region length":True if self.targetsInRegions or self.gffFile else False, "Alignment strand":True, "Start block position in target":True, "End block position in target":True, "Number of sirnas in block":True, "Total sirnas in target region" if self.targetsInRegions else "Total sirnas in target":True, "Total blocks in target region" if self.targetsInRegions else "Total blocks in target":True, "Target region annotation" if self.targetsInRegions else "Target annotation":True if self.gffFile else False}
            hoverName = "Target ID"

        fig = px.timeline(df, x_start="Start position in query", x_end="End position in query", y="Target region name", color="Target ID", hover_name=hoverName, hover_data=hoverData)
        fig.layout.xaxis.type = 'linear'
        #Modify x-labels to linear
        for figData in fig.data:
            xData = []
            yNamePre = None
            pos = 0
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
        return self.jsonFileName
    
    def dataToJson(self):
        #Extracts the data from bowtie results and efficiency and put everything into json format.
        for sirnaName in self.sirnaData:
            sirna = self.sirnaData[sirnaName]
            if sirna.bowtieDataTuple():
                for alignment in sirna.bowtieDataTuple():
                    hitTarget = alignment[0]
                    offTarget = hitTarget not in self.mainTargets
                    sirnaDict = {"query_name": self.queryName, "sirna_position": sirnaName}
                    sirnaDict.update(sirna.toDict())
                    sirnaDict.update({"is_off_target": offTarget,
                                      "hit_name": hitTarget.getName(),
                                      "hit_region":"-".join(map(str,hitTarget.getRegion())),
                                      "hit_annotation":" | ".join(hitTarget.getAnnotation()),
                                      "reference_pos": alignment[1],
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

