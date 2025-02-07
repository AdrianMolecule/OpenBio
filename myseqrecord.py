from re import S
from Bio.Seq import Seq
from typing import Optional
from typing import TypeVar
from typing import Union

from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from Bio.SeqFeature import SeqFeature, FeatureLocation
# Define a TypeVar for the class type
MyT = TypeVar('T', bound='MySeqRecord')

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
		self.notAnealedLocation:tuple[int, int]=None
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
		collectDict['notAnealedLocation'] =  str(self.notAnealedLocation)
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

	def setNotAnealedLocation(self, startStop:tuple):
		# if not self.isPrimer:
		# 	raise ValueError(f"Invalid operation:  The aneal location can be set only for the primers not for the strands")
		self.notAnealedLocation: tuple[int, int]=startStop
		

	def getNotAnealedLocation(self)->tuple:
		return self.notAnealedLocation
		
	def __repr__(self):
		# Customize the representation of the item
		return "id"+str(self.id)+" description:"+self.description+" SingleStranded:"+str(self.singleStranded)+" primer:"+str(self.isPrimer)+" 5To3:"+str(self.fiveTo3) +" primer:"+str(self.isPrimer)+" Sequence:"+str( self.seq)[:70]


	def addFeature(self, start, end, strand:int, type, id, label:str):
		newFeature:SeqFeature=SeqFeature(FeatureLocation(start, end, strand=None), type=type, id=id, qualifiers={"label": [label]  })
		self.features.append(newFeature)   


	def removeFullSpanningFeatures(self): # that is usually where type=="source"
		for i, feature in enumerate(self.features):	
			if self.toIgnore(feature):
				self.features.pop(i)		


	def toIgnore(self,feature:SeqFeature)->bool: # usually for record.type=="source" that takes the full length of the sequence):
		if feature.location.start==0 and feature.location.end==len(self):
			return True
			return False

	def splitRecord(self, splitPointIndex, extraIndentForSecond=0)->MyT:
		if splitPointIndex>=len(self.seq):
			return
		seqString=MySeqRecord.seqToString(self.seq)
	
		seqU:Seq=  Seq(seqString[splitPointIndex:])
		seqL:Seq=  Seq(seqString[0:splitPointIndex])
		newUpper=MySeqRecord(SeqRecord(seqU,id=self.id+"_t", name=self.name, description="truncated"+ self.description), singleStranded=True,fiveTo3=True,primer=False)
		newUpperFeatures=self.shiftFeaturesLocs(splitPointIndex)
		newUpper.features=newUpperFeatures
		newUpper.xStartOffsetAsLetters=self.xStartOffsetAsLetters+splitPointIndex
		newLower=MySeqRecord(SeqRecord(seqL,id=self.id+"_l", name=self.name, description="looped"+ self.description), singleStranded=True,fiveTo3=True,primer=True)

		newLower.xStartOffsetAsLetters=newUpper.xStartOffsetAsLetters+extraIndentForSecond
		return newUpper, newLower
	
		# for a in self.annotations:
		# 	start=a.location.start
		# 	end=a.location.end
		# 	if start>byHowMuch:
		# 		if end >byHowMuch:
		# 			a.location.start-=byHowMuch
		# 			a.location.end-=byHowMuch
		# 		else:
		# 			a.location.start=byHowMuch
		# 			a.location.end-=byHowMuch


	def shiftFeaturesLocs(self, byHowMuch)->list[SeqFeature]:
		newFeatures=list()
		for f in self.features:
			f:SeqFeature			
			start=f.location.start
			end=f.location.end
			if start>byHowMuch:				
				if end >byHowMuch:
					l=FeatureLocation(f.location.start-byHowMuch,f.location.end-byHowMuch, strand=+1)
				else:
					l=FeatureLocation(byHowMuch,f.location.end-byHowMuch, strand=+1)
				newFeature:SeqFeature=SeqFeature(l, type=f.type, id=f.id, qualifiers={"label": ["High"+f.qualifiers.get("label")[0]]  })
				newFeatures.append(newFeature)
		return newFeatures
	
	def shiftBackwardsFeaturesLocs(self, lowerLoopLen)->list[SeqFeature]:
		newFeatures=list()
		for f in reversed(self.features):
			f:SeqFeature	

			end=lowerLoopLen-f.location.start
			start=lowerLoopLen-f.location.end
			if start>0:				
				if start >0:# fits the hole feature 
					l=FeatureLocation(start,end, strand=-1)
				# else:
				# 	l=FeatureLocation(byHowMuch,f.location.end-byHowMuch, strand=-1)
				newFeature:SeqFeature=SeqFeature(l, type=f.type, id=f.id, qualifiers={"label": ["Lo"+f.qualifiers.get("label")[0]]  })
				newFeatures.append(newFeature)
		return newFeatures	
			
	def shiftLocations(self, byHowMuch):
		newFeatures=list()
		for f in self.features:
			f:SeqFeature			
			start=f.location.start
			end=f.location.end
			if start>byHowMuch:				
				if end >byHowMuch:
					l=FeatureLocation(f.location.start-byHowMuch,f.location.end-byHowMuch, strand=+1)
				else:
					l=FeatureLocation(byHowMuch,f.location.end-byHowMuch, strand=+1)
				newFeature:SeqFeature=SeqFeature(l, type=f.type, id=f.id, qualifiers={"label": ["YY"+f.qualifiers.get("label")[0]]  })
				newFeatures.append(newFeature)
		self.features=newFeatures
		# for a in self.annotations:
		# 	start=a.location.start
		# 	end=a.location.end
		# 	if start>byHowMuch:
		# 		if end >byHowMuch:
		# 			a.location.start-=byHowMuch
		# 			a.location.end-=byHowMuch
		# 		else:
		# 			a.location.start=byHowMuch
		# 			a.location.end-=byHowMuch	

	def seqToString(seq:Seq)->str:
		return seq._data.decode('ASCII') # from bytes to stringead and strip any extra whitespace or newline characters
