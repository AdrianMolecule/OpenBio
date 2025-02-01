from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord
from myseqrecord import MySeqRecord

class PrimerUtils:
    # return a list of tuples each containig a start and an end
    def findPrimerOverlaps(targetDnaRecordSequence:Seq, primerRecordSequence:Seq, minOverlapLength):
        overlapStartsEnds = set()
        largestOverlaps = []
        largestLength = 0
        largestOverlapsInShorts = []
        # if targetDnaRecord.is
        from util import seqToString
        complementdPrimerString:str=seqToString(primerRecordSequence).lower()
        targetDnaRecordString: str =seqToString(targetDnaRecordSequence).lower()

        for i in range(len(complementdPrimerString) - minOverlapLength + 1):
            for j in range(i + minOverlapLength, len(complementdPrimerString) + 1):
                subShort: str = complementdPrimerString[i:j]
                index: int = targetDnaRecordString.find(subShort)
                while index != -1:
                    overlapStartsEnds.add((index, index + len(subShort) - 1, i, j - 1))
                    if len(subShort) > largestLength:
                        largestLength = len(subShort)
                        largestOverlaps: list[tuple[int, int]] = [(index, index + len(subShort) - 1)]
                        largestOverlapsInShorts: list[tuple[int, int]] = [(i, j - 1)]
                    elif len(subShort) == largestLength:
                        largestOverlaps.append((index, index + len(subShort) - 1))
                        largestOverlapsInShorts.append((i, j - 1))
                    index = targetDnaRecordString.find(subShort, index + 1)

        return list(overlapStartsEnds), largestOverlaps, largestOverlapsInShorts


def main():
    targetDnaSequence:Seq =         Seq("aaAATTccGGCCa")# 3To5
    uncomplementedPrimerRecord:Seq =Seq("ttttAATTaGGCt") # and this one complemented looks like aaaaTTAAtCCGa
    complementedPrimerSeq=uncomplementedPrimerRecord.complement()

    minOverlapLength = 3
    overlaps, largestInStrand, largestInPrimer = PrimerUtils.findPrimerOverlaps(targetDnaSequence, complementedPrimerSeq, minOverlapLength)
    print("targetDnaSequence:", targetDnaSequence)
    print("uncomplementedPrimerRecord:", uncomplementedPrimerRecord)
    print("primerSeq:", complementedPrimerSeq)
    print("Overlaps:", overlaps)
    print("largest:", largestInStrand)
    from util import seqToString
    if largestInStrand:
        print("Largest Overlaps in targetDnaSequence:", largestInStrand, seqToString(targetDnaSequence)[largestInStrand[0][0]:largestInStrand[0][1]] )
        print("Largest Overlaps in primer:", largestInPrimer, seqToString(targetDnaSequence)[largestInPrimer[0][0]:largestInPrimer[0][1]] )


# if __name__ == "__main__":
#     main()


    
