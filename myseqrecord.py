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
        collectDict = {}            
        # collectDict the attributes from the parent SeqRecord class
        collectDict['id'] = str(self.id)
        collectDict['name'] = str(self.name)
        collectDict['description'] = self.description
        
        # collectDict additional attributes from MySeqRecord
        collectDict['singleStranded'] = str(self.singleStranded)
        collectDict['isPrimer'] =  str(self.isPrimer)
        collectDict['fiveTo3'] =  str(self.fiveTo3)
        collectDict['xStartOffsetAsLetters'] =  str(self.xStartOffsetAsLetters)
        if self.hybridizedToStrand:
            collectDict['hybridizedToStrand'] =  self.hybridizedToStrand.description+"-"+str("5 to 3 " if self.fiveTo3 else "3 to 5")
        if self.hybridizedToPrimer:
            collectDict['hybridizedToPrimer'] =  self.hybridizedToPrimer.description
        collectDict['uniqueId'] =  str(self.uniqueId)
        collectDict['seq'] = self.seq._data.decode('ASCII')
        featuresString="\n"
        for f in self.features:
            featuresString+=f.type+"   "
            if f.qualifiers.get("label"):
                featuresString+=str(f.location)+f.qualifiers.get("label")[0]
                featuresString+="\n"
        collectDict['features'] = featuresString
        annotationsString="\n"
        #type id qualifier
        for key, value in self.annotations.items():
            annotationsString+=(f"{key}= {value}")        
            annotationsString+="\n"
        collectDict['annotations'] = annotationsString        
        collectDict['dbxrefs'] =  str(self.dbxrefs)
        s=""
        for key in collectDict:
            s+=key+":"+collectDict[key]+"\n"
        return s

    def __repr__(self):
        # Customize the representation of the item
        return "id"+str(self.id)+" description:"+self.description+" SingleStranded:"+str(self.singleStranded)+" primer:"+str(self.isPrimer)+" 5To3:"+str(self.fiveTo3) +" primer:"+str(self.isPrimer)+" Sequence:"+str( self.seq)[:70]


        
