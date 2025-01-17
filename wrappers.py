from re import S
from Bio.Seq import Seq
from typing import Optional

from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from typing import Union

class MySeqRecord(SeqRecord):

    def __init__(self, seqRecord):
        # Call the parent class constructor with the attributes of seqRecord
        super().__init__(seqRecord.seq, id=seqRecord.id, name=seqRecord.name, description=seqRecord.description,  annotations=seqRecord.annotations, features=seqRecord.features, dbxrefs=seqRecord.dbxrefs)
        # Copy all attributes from the original SeqRecord
        for key, value in vars(seqRecord).items():
            setattr(self, key, value)  
            # 2 added fields      
        self.singleStranded = False
        self.shrinkedFeatures=None
        self.isPrimer=False
        

    def displayInfo(self):
        print(f"ID: {self.id}")
        print(f"Description: {self.description}")
        print(f"Sequence: {self.seq}")
        print(f"Extra Info: {self.extra_info}")

        
