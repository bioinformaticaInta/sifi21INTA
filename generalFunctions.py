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
