from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord

from util import *
from wrappers import MySeqRecord

def findPrimerOverlaps(targetDnaRecord:MySeqRecord, primerRecord:MySeqRecord, minOverlapLength):
    overlapStartsEnds = set()
    largestOverlaps = []
    largestLength = 0
    largestOverlapInShorts = []
    # if targetDnaRecord.is
    complementdPrimerString:str=seqToString(primerRecord.seq.complement()).lower()
    targetDnaRecordString: str =seqToString(targetDnaRecord.seq).lower()

    for i in range(len(complementdPrimerString) - minOverlapLength + 1):
        for j in range(i + minOverlapLength, len(complementdPrimerString) + 1):
            subShort: str = complementdPrimerString[i:j]
            index: int = targetDnaRecordString.find(subShort)
            while index != -1:
                overlapStartsEnds.add((index, index + len(subShort) - 1, i, j - 1))
                if len(subShort) > largestLength:
                    largestLength = len(subShort)
                    largestOverlaps: list[tuple[int, int]] = [(index, index + len(subShort) - 1)]
                    largestOverlapInShorts: list[tuple[int, int]] = [(i, j - 1)]
                elif len(subShort) == largestLength:
                    largestOverlaps.append((index, index + len(subShort) - 1))
                    largestOverlapInShorts.append((i, j - 1))
                index = targetDnaRecordString.find(subShort, index + 1)

    return list(overlapStartsEnds), largestOverlaps, largestOverlapInShorts

longString =  MySeqRecord(SeqRecord("aaAATTccGGCCa"),True,fiveTo3=False,primer=False)
shortString =  MySeqRecord(SeqRecord("aaaaTTAAtCCGa"),True,fiveTo3=True,primer=True)
minOverlapLength = 3
overlaps, largest, largestInShort = findPrimerOverlaps(longString, shortString, minOverlapLength)
print("Overlaps:", overlaps)
print("Largest Overlaps in longString:", largest)
print("Largest Overlaps in shortString:", largestInShort)


    
