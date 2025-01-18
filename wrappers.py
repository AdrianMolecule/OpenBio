from re import S
from Bio.Seq import Seq
from typing import Optional

from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from typing import Union

class MySeqRecord(SeqRecord):

    def __init__(self, seqRecord,singleStranded:bool,fiveTo3:bool,primer:bool):
        # Call the parent class constructor with the attributes of seqRecord
        super().__init__(seqRecord.seq, id=seqRecord.id, name=seqRecord.name, description=seqRecord.description,  annotations=seqRecord.annotations, features=seqRecord.features, dbxrefs=seqRecord.dbxrefs)
        # Copy all attributes from the original SeqRecord
        for key, value in vars(seqRecord).items():
            setattr(self, key, value)  
            # 2 added fields      
        self.singleStranded = False
        self.fiveTo3 = True
        self.primer=False
        self.hybridizedTo:MySeqRecord=None
        self.shrinkedFeatures=None

    def __str__(self):
            # Customize how the item is printed
            return f"Item: {self.value}"

    def __repr__(self):
        # Customize the representation of the item
        return "id"+str(self.id)+" description:"+self.description+" SingleStranded:"+str(self.singleStranded)+ " primer:"+str(self.primer)+" 5To3:"+self.fiveTo3+" primer:"+self.primer,+" Sequence:"+str( self.seq)

        

    def displayInfo(self):
        print(f"ID: {self.id}")
        print(f"Description: {self.description}")
        print(f"Sequence: {self.seq}")
        print(f"Extra Info: {self.extra_info}")

        
