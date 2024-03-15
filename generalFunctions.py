#!/usr/bin/python

from BCBio import GFF
from datetime import datetime
import sys
import sifi21INTA.TargetSequence as TargetSequence

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

def addTargetsRegions(bowtieAlignments, gapSize, xmer):
    #The input is a list of lists containing: [sirna, hitTarget, hitPos, hitStrand, hitMissmatches]
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
                    #Change the end of the region for the new end (extend) and set the flag
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
    gffData = list(GFF.parse(gffFile))
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
        annotation = targetsAnnotations[hitTarget.getName()]
        if annotation:
            hitTarget.addAnnotation(annotation[0])
            hitTarget.setLength(annotation[1])
    
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
    refInfo = tuple()
    for rec in gffData:
        for feature in rec.features:
            for subfeature in feature.sub_features:
                if seqName in subfeature.id or ('Name' in subfeature.qualifiers and seqName in subfeature.qualifiers['Name']):
                    refInfo = (feature.qualifiers["description"][0],subfeature.location.end-subfeature.location.start+1)
    return refInfo
