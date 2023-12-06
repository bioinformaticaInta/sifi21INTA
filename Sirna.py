#!/usr/bin/python

from Bio.Seq import Seq

import sifi21INTA.freeEnergy as freeEnergy

class Sirna:
    def __init__(self, sequence):
        self.sequence = sequence
        self.bowtieData = []     #list of tuples=(hitName, hitPos, hitStrand, hitMissmatches)
        self.rnaplfoldData = tuple()
        self.lunpDataXmer = None
        self.strandSelection = None
        self.endStability = None
        self.targetSiteAccessibility = None
        self.isThermoEfficient = None
        self.isEfficient = None
   
    def getSequence(self):
        return self.sequence

    def bowtieDataTuple(self):
        return self.bowtieData

    def __repr__(self):
        return '< ' + ' | '.join(["sequence="+self.sequence,"bowtieData="+str(self.bowtieData)]) + ' >'

    def toDict(self):
        return {"sirna_sequence": self.sequence, 
                "is_efficient": self.isEfficient,
                "is_thermo_efficient": self.isThermoEfficient,
                "strand_selection": self.strandSelection, 
                "end_stability": self.endStability,
                "target_site_accessibility": self.targetSiteAccessibility,
                "accessibility_value": self.lunpDataXmer, 
                }
    
    def getEfficiency(self):
        return self.isEfficient == True

    def bowtieHitNames(self):
        return {alignment[0] for alignment in self.bowtieData} if self.bowtieData else set() 

    def getOffTargets(self, mainTargets):
        offTargets = set()
        for hitName in self.bowtieHitNames():
            if hitName not in mainTargets:
                offTargets.add(hitName)
        return offTargets

    def addBowtieAlignment(self, name, pos, strand, missmatches):
        self.bowtieData.append((name, pos, strand, missmatches)) 

    def addRnaplfoldData(self, lunpData):
        self.rnaplfoldData = tuple(map(float, lunpData))    

    def calculateEfficiency(self, sequenceN2, accessibilityWindow, tsAccessibilityTreshold, endStabilityTreshold, startPosition, endNucleotides, overhang, terminalCheck, strandCheck, endCheck, accessibilityCheck):
        self.lunpDataXmer = self.rnaplfoldData[accessibilityWindow-1]
        self.calculateStrandSelection(sequenceN2, startPosition, endNucleotides, overhang)
        self.calculateEndStability(sequenceN2, endStabilityTreshold, startPosition, endNucleotides, overhang)
        self.calculatePairProbability(tsAccessibilityTreshold)
        self.checkThermoEfficient(strandCheck, endCheck, accessibilityCheck)

        sirnaSize = len(self.sequence)
        if terminalCheck:
            if self.sequence[sirnaSize-3] == 'A' or self.sequence[sirnaSize-3] == 'T':
                #print "A/T at S19"
                if self.sequence[1] == 'A' or self.sequence[1] == 'T':
                    #print "A/T at S1"
                    if self.isThermoEfficient:
                        self.isEfficient = True
                    else:
                        self.isEfficient = False
                else:
                    if accessibilityCheck:
                        if self.targetSiteAccessibility:
                            self.isEfficient = True
                        else:
                            self.isEfficient = False
                    else:
                        self.isEfficient = True
            else:
                if self.sequence[1] == 'G' or self.sequence[1] == 'C':
                    #print "C/G at S1"
                    if self.isThermoEfficient:
                        self.isEfficient = True
                    else:
                        self.isEfficient = False
                else:
                    self.isEfficient = False
        else:
            if self.isThermoEfficient:
                self.isEfficient = True
            else:
                self.isEfficient = False

    def calculateStrandSelection(self, sequenceN2, startPosition, endNucleotides, overhang):
        """Returns whether the strand will be selected (True) or not (False) based on energy rules."""
        # For siRNA n>3 we calculate with dangling ends
        if sequenceN2 != None:
            sense5_MFE_enegery, anti_sense5_MFE_enegery = self.freeEnergyDanglingEnds(sequenceN2, startPosition, endNucleotides)
        else:
            # For the first two siRNAs no dangling ends
            sense5_MFE_enegery, anti_sense5_MFE_enegery = self.freeEnergy(startPosition, endNucleotides, overhang)

        if anti_sense5_MFE_enegery >= sense5_MFE_enegery:
            self.strandSelection = True
        else:
            self.strandSelection = False
        #print "G_S/G_AS ", sense5_MFE_enegery, anti_sense5_MFE_enegery

    def calculateEndStability(self, sequenceN2, endStabilityTreshold, startPosition, endNucleotides, overhang):
        """Calculate whether the end stability is higher or equal threshold (default=1).
           Return True if yes and False if it is lower."""
        # For siRNA n>3 we calculate with dangling ends
        if sequenceN2 != None:
            sense5_MFE_enegery, anti_sense5_MFE_enegery = self.freeEnergyDanglingEnds(sequenceN2, startPosition, endNucleotides)
        else:
            # For the first two siRNAs no dangling ends
            sense5_MFE_enegery, anti_sense5_MFE_enegery = self.freeEnergy(startPosition, endNucleotides, overhang)

        # End stability
        if (anti_sense5_MFE_enegery - sense5_MFE_enegery) >= endStabilityTreshold:
            self.endStability = True
        else:
            self.endStability = False

    def calculatePairProbability(self, tsAccessibilityTreshold):
        """Calculates whether the pair probability the siRNA at a certain window (default 8) is higher or equal
           the threshold (default=0.1)
           Return True if yes and False if it is lower."""
        # Get pair probability of the accessibility window (chosen by user) of siRNA
        if self.lunpDataXmer >= tsAccessibilityTreshold:
            self.targetSiteAccessibility = True
        else:
            self.targetSiteAccessibility = False

    def checkThermoEfficient(self, strandCheck, endCheck, accessibilityCheck):
        """Calculates whether a siRNA is efficient or not. Default priority rules (if checked by user):
           1. Strand selection must be True.
           2. End stability must be higher or equal threshold (default=1).
           3. Pair probability (lunp data) of xmer (chosen by user, default 8) must be higher or equal threshold (default=0.1).
           If all rules apply, a siRNA is efficient."""

        if strandCheck and endCheck and accessibilityCheck:
            if self.strandSelection and self.endStability and self.targetSiteAccessibility:
                self.isThermoEfficient = True
            else:
                self.isThermoEfficient = False
        if strandCheck and endCheck and not accessibilityCheck:
            if self.strandSelection and self.endStability:
                self.isThermoEfficient = True
            else:
                self.isThermoEfficient = False
        if strandCheck and accessibilityCheck and not endCheck:
            if self.strandSelection and self.targetSiteAccessibility:
                self.isThermoEfficient = True
            else:
                self.isThermoEfficient = False
        if endCheck and accessibilityCheck and not strandCheck:
            if self.endStability and self.targetSiteAccessibility:
                self.isThermoEfficient = True
            else:
                self.isThermoEfficient = False
        if strandCheck and not endCheck and not accessibilityCheck:
            if self.strandSelection:
                self.isThermoEfficient = True
            else:
                self.isThermoEfficient = False
        if endCheck and not strandCheck and not accessibilityCheck:
            if self.endStability:
                self.isThermoEfficient = True
            else:
                self.isThermoEfficient = False
        if accessibilityCheck and not strandCheck and not endCheck:
            if self.targetSiteAccessibility:
                self.isThermoEfficient = True
            else:
                self.isThermoEfficient = False
        if not accessibilityCheck and not strandCheck and not endCheck:
            self.isThermoEfficient = True

    def freeEnergy(self, startPosition, endNucleotides, overhang):
        """Calculate the free energy of a sequence.

           Code was taken from Biopython.
           http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-pysrc.html

           Example for 21mer
           siRNA GGGATGGCTCAAAGGCGTAGT
           Sense5prime_MFE siRNA position [0,1,2] -> GGG
           Antisense5prime_MFE siRNA position [17,18,29] -> GTA"""

        sirnaSize = len(self.sequence)
        # Sense5_MFE
        sense_five_seq = self.sequence[startPosition : endNucleotides]
        # Anitsense5_MFE
        antisense_five_seq = self.sequence[sirnaSize-overhang-endNucleotides : sirnaSize-overhang]

        sense5_MFE_enegery = freeEnergy.calculateFreeEnergy(sense_five_seq)
        anti_sense5_MFE_enegery = freeEnergy.calculateFreeEnergy(antisense_five_seq)

        return sense5_MFE_enegery, anti_sense5_MFE_enegery

    def freeEnergyDanglingEnds(self, sequenceN2, startPosition, endNucleotides):
        """Calculate the free energy of a sequence.
           c_seq is for dangling ends"""

        sirnaSize = len(self.sequence)
        #Sense5_MFE
        sense_five_seq = self.sequence[startPosition : endNucleotides]
        sense_c_seq = Seq(sequenceN2).reverse_complement().strip()[sirnaSize-5 : sirnaSize-1]
        # Anitsense5_MFE
        antisense_five_seq = Seq(sequenceN2).reverse_complement().strip()[startPosition : endNucleotides]
        antisense_c_seq = self.sequence[sirnaSize-5 : sirnaSize-1]

        sense5_MFE_enegery = freeEnergy.calculateFreeEnergy(sense_five_seq, check=True, strict=True, c_seq=sense_c_seq[::-1], shift=1)
        anti_sense5_MFE_enegery = freeEnergy.calculateFreeEnergy(antisense_five_seq, check=True, strict=True, c_seq=antisense_c_seq[::-1], shift=1)
        #print 'G ', anti_sense5_MFE_enegery

        return sense5_MFE_enegery, anti_sense5_MFE_enegery
