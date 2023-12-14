#!/usr/bin/python

import re
import json
from collections import Counter

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

def addChromosomeRegions(bowtieAlignments, gapSize, xmer):
    #The input is a list of lists containing: [sirna, hitName, hisPos, hitStrand, hitMissmatches]
    chromRegions = {}
    #Traversing ordered by hit positions assembling the regions shared between matches
    for alignment in sorted(bowtieAlignments, key =lambda x: x[2]):
        start = alignment[2]
        end = start + xmer-1
        if not alignment[1] in chromRegions:
            chromRegions[alignment[1]] = [[start,end]]
        else:
            #if any region exists, search the possible overlap
            regions = chromRegions[alignment[1]]
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
        regions = chromRegions[alignment[1]]
        for region in regions:
            if region[0] <= alignment[2] <= region[1]:
                alignment[1] += "_"+"-".join(list(map(str,region)))
            
def bowtieToList(bowtieFileName):
    bowtieFile = open(bowtieFileName)
    bowtieAlignments = []                                                                                                              
    for bowtieMatch in bowtieFile:
        bowtieMatch = bowtieMatch.strip().split('\t')                                                                                  
        sirnaName = int(bowtieMatch[0])
        hitStrand = bowtieMatch[1]
        hitName = bowtieMatch[2]
        hitPos = int(bowtieMatch[3])+1
        hitMissmatches = 0
        if len(bowtieMatch) == 8:
            hitMissmatches = int(bowtieMatch[7])
        bowtieAlignments.append([sirnaName, hitName, hitPos, hitStrand, hitMissmatches])
    bowtieFile.close()
    return bowtieAlignments
