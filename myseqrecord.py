from re import S
from Bio.Seq import Seq
from typing import Optional
from typing import Union

from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from Bio.SeqFeature import SeqFeature, FeatureLocation

class MySeqRecord(SeqRecord):
    uniqueId=0
    def __init__(self, seqRecord,singleStranded:bool,fiveTo3:bool,primer:bool):
        # Call the parent class constructor with the attributes of seqRecord
        super().__init__(seqRecord.seq, id=seqRecord.id, name=seqRecord.name, description=seqRecord.description,  annotations=seqRecord.annotations, features=seqRecord.features, dbxrefs=seqRecord.dbxrefs)
        # Copy all attributes from the original SeqRecord
        for key, value in vars(seqRecord).items():
            setattr(self, key, value)  
        self.singleStranded = singleStranded
        self.isPrimer=primer
        self.fiveTo3 = fiveTo3
        self.xStartOffsetAsLetters=0 # in letters not pixels       
        self.hybridizedToStrand:MySeqRecord=None
        self.hybridizedToPrimer:MySeqRecord=None
        self.uniqueId=MySeqRecord.uniqueId
        MySeqRecord.uniqueId+=1

    def __str__(self):
            # Customize how the item is printed
            return f"Item: {self.value}"

    def __repr__(self):
        # Customize the representation of the item
        return "id"+str(self.id)+" description:"+self.description+" SingleStranded:"+str(self.singleStranded)+" primer:"+str(self.isPrimer)+" 5To3:"+str(self.fiveTo3) +" primer:"+str(self.isPrimer)+" Sequence:"+str( self.seq)[:70]

    def displayInfo(self):
        print(f"ID: {self.id}")
        print(f"Description: {self.description}")
        print(f"Sequence: {self.seq}")
        print(f"Extra Info: {self.extra_info}")

        
