#!/usr/bin/python

import re
import json
from collections import Counter
from BCBio import GFF

import sifi21INTA.TargetSequence as TargetSequence

def iterparse(j):
    """Work around because json can not load multiple objects.
       Taken from http://stackoverflow.com/questions/22112439/valueerror-extra-data-while-loading-json."""
    nonspace = re.compile(r'\S')
    decoder = json.JSONDecoder()
    pos = 0
    while True:
        matched = nonspace.search(j, pos)
        if not matched:
            break
        pos = matched.start()
        decoded, pos = decoder.raw_decode(j, pos)
        yield decoded


def prepareJsonData(jsonFileName):
    """Prepares the json file for easier parsing."""
    sifi_data = open(jsonFileName, "r").read()
    data = list(iterparse(sifi_data))[0]
    return data

def getTableData(jsonFileName):
    """Extracts a summary of all and efficient hits."""
    query = prepareJsonData(jsonFileName)
    # Get number of hits per query
    hitCounter = Counter(player['hit_name'] for player in query)
    efficicentCounter = Counter(player['hit_name'] for player in query if player['is_efficient'])
    tableData = []
    for x in hitCounter.most_common():
        for y in efficicentCounter.most_common():
            if x[0] == y[0]:
                tableData.append([x[0], x[1], y[1]])
        if x[0] not in list(efficicentCounter):
            tableData.append([x[0], x[1], 0])
    return tableData

def getChromosomeRegion(pos, chromLen, regionSize):
    regionStr = None
    if 1 <= pos <= chromLen:
        regionStr = "%s-%s" % searchInterval(tuple(range(1,chromLen+1, regionSize)), pos)
    return regionStr

#Binary search of positional range interval
def searchInterval(tuplePos, pos):
    region = (None,None)
    tupleLen = len(tuplePos)
    if tupleLen > 2:
        center = len(tuplePos)//2
        if pos < tuplePos[center]:
            if pos >= tuplePos[center-1]:
                region = (tuplePos[center-1], tuplePos[center])
            else:
                region = searchInterval(tuplePos[:center], pos)
        else:
            if pos < tuplePos[center+1]:
                region = (tuplePos[center], tuplePos[center+1])
            else:
                region = searchInterval(tuplePos[center+1:], pos)
    return region

def addTargetsRegions(bowtieAlignments, gapSize, xmer):
    #The input is a list of lists containing: [sirna, hitTarget, hisPos, hitStrand, hitMissmatches]
    targetsRegions = {}
    #Traversing ordered by hit positions assembling the regions shared between matches
    for alignment in sorted(bowtieAlignments, key =lambda x: x[2]):
        start = alignment[2]
        end = start + xmer-1
        hitTarget = alignment[1]
        if not hitTarget.getName() in targetsRegions:
            targetsRegions[hitTarget.getName()] = [[start,end]]
        else:
            #if any region exists, search the possible overlap
            regions = targetsRegions[hitTarget.getName()]
            posRegion = 0
            prevRegion = False
            while not prevRegion and posRegion < len(regions):
                region = regions[posRegion]
                if region[0] <= start <= region[1] or (start-region[1]) <= gapSize: #overlap!!!
                    region[1] = end
                    prevRegion = True
                posRegion += 1
            if not prevRegion: #Insert new regions if not overlaping
                regions.append([start,end])
    #Add position of the region to each hit name
    for alignment in bowtieAlignments:
        hitTarget = alignment[1]
        regions = targetsRegions[hitTarget.getName()]
        for region in regions:
            if region[0] <= alignment[2] <= region[1]:
                hitTarget.setRegion(*region)
    return targetsRegions

def bowtieToList(bowtieFileName):
    bowtieFile = open(bowtieFileName)
    bowtieAlignments = []                                                                                                              
    targetsRegions = {}
    for bowtieMatch in bowtieFile:
        bowtieMatch = bowtieMatch.strip().split('\t')                                                                                  
        sirnaName = int(bowtieMatch[0])
        hitStrand = bowtieMatch[1]
        hitTarget = TargetSequence.TargetSequence(bowtieMatch[2])
        hitPos = int(bowtieMatch[3])+1
        hitMissmatches = 0
        if len(bowtieMatch) == 8:
            hitMissmatches = int(bowtieMatch[7])
        bowtieAlignments.append([sirnaName, hitTarget, hitPos, hitStrand, hitMissmatches])
        targetsRegions[bowtieMatch[2]] = []
    bowtieFile.close()
    return bowtieAlignments,targetsRegions

def addAnnotations(gffFileName, bowtieAlignments, genomicRef, inRegions, targetsRegions):
    if genomicRef and inRegions:
        gffData = getGFFDataFromFile(gffFileName)
        #Add annotations from GFF to each region
        #targetRegions is a dictionary with {targetName:[targetRegion1, targetRegion2, etc]} / targetRegion1=[startRegion1,endRegion1]
        addGenomicAnnotations(gffData, bowtieAlignments, targetsRegions)
    elif not genomicRef:
        gffData = getGFFDataFromFile(gffFileName)
        #targetsRegions is a dictionary as the genomic case with regions as value (targetsInRegions True) or with empty list (targetsInRegions False)
        addTranscriptomicAnnotations(gffData, bowtieAlignments, targetsRegions.keys())

def getGFFDataFromFile(gffFileName):
    gffFile = open(gffFileName)
    gffData = []
    for rec in GFF.parse(gffFile):
        gffData.append(rec)
    gffFile.close()
    return gffData 

def addGenomicAnnotations(gffData, bowtieAlignments, targetsRegions):
    targetsAnnotations = {}
    for targetName in targetsRegions:
        for region in targetsRegions[targetName]:
            targetsAnnotations[(targetName, *region)] = getGenomicRegionAnnotation(gffData, targetName, *region)
    #Add annotation for each alignment
    for alignment in bowtieAlignments:
        hitTarget = alignment[1]
        for annotation in targetsAnnotations[(hitTarget.getName(), *hitTarget.getRegion())]:
            hitTarget.addAnnotation(annotation)

def addTranscriptomicAnnotations(gffData, bowtieAlignments, targetsNames):
    targetsAnnotations = {}
    for targetName in targetsNames:
        targetsAnnotations[targetName] = getTranscriptAnnotation(gffData, targetName)
    #Add annotation for each alignment
    for alignment in bowtieAlignments:
        hitTarget = alignment[1]
        for annotation in targetsAnnotations[hitTarget.getName()]:
                hitTarget.addAnnotation(annotation)
    
def getGenomicRegionAnnotation(gffData, targetName, startRegion, endRegion):
    refInfo = []
    for rec in gffData:
        if rec.id == targetName:
            for feature in rec.features: ##agregar condicion para que solo sean genes
                if feature.type == "gene":
                    featureInRegion = None
                    #Partial end of gene at the region start
                    if feature.location.start < startRegion and startRegion < feature.location.end <= endRegion:
                        featureInRegion = feature
                    #Partial start of gene at the region end
                    elif startRegion <= feature.location.start < endRegion and  feature.location.end > endRegion:
                        featureInRegion = feature
                    #gene included in region
                    elif feature.location.start >= startRegion and  feature.location.end <= endRegion:
                        featureInRegion = feature
                    #region included in gene
                    elif feature.location.start < startRegion and  feature.location.end > endRegion:
                        featureInRegion = feature
                    if featureInRegion:
                        data = featureInRegion.id
                        if 'description' in featureInRegion.qualifiers:
                            data += " (" + featureInRegion.qualifiers['description'][0] + ")"
                        refInfo.append(data)
    return refInfo

def getTranscriptAnnotation(gffData, seqName):
    refInfo = []
    for rec in gffData:
        for feature in rec.features:
            for subfeature in feature.sub_features:
                if seqName in subfeature.id or ('Name' in subfeature.qualifiers and seqName in subfeature.qualifiers['Name']):
                    refInfo.append(feature.qualifiers["description"][0])
    return refInfo

