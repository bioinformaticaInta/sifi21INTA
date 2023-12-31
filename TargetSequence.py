#!/usr/bin/python

class TargetSequence:
    def __init__(self, name, startRegion=None, endRegion=None, length=None, annotation=None):
        self.name = name
        self.startRegion = startRegion
        self.endRegion = endRegion
        self.length = length
        self.annotation = []
   
    def __repr__(self):
        return '< ' + ' | '.join(["ID="+self.name,"Region="+str(self.startRegion)+"-"+str(self.endRegion),"Length="+str(self.length),"Annotation="+str(self.annotation)]) + ' >'
    
    def __eq__(ts1, ts2):
        return ts1.name==ts2.name and ts1.startRegion==ts2.startRegion and ts1.endRegion==ts2.endRegion

    def __hash__(self):
        return(hash(self.name+str(self.startRegion)+str(self.startRegion)))

    def getName(self):
        return self.name

    def getRegion(self):
        return self.startRegion,self.endRegion
    
    def getNameWithRegion(self):
        name = self.name
        if self.startRegion and self.endRegion:
            name += "_"+str(self.startRegion)+"-"+str(self.endRegion)
        return name
    
    def getAnnotation(self):
        return self.annotation

    def getLength(self):
        return self.length

    def setRegion(self, startRegion, endRegion):
        self.startRegion = startRegion
        self.endRegion = endRegion
        self.length = endRegion-startRegion+1

    def setLength(self, length):
        if not self.startRegion and not self.endRegion:
            self.length = length

    def addAnnotation(self, annotation):
        self.annotation.append(annotation)

    def isRegion(self, start, end):
        return start==self.startRegion and end==self.endRegion
