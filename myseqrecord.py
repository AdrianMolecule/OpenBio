from re import S
from Bio.Seq import Seq
from typing import Optional
from typing import TypeVar
from typing import Union

from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from Bio.SeqFeature import SeqFeature, SimpleLocation


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
		self.notAnnealedLocation:tuple[int, int]=None
		self.loopInfo:tuple[MySeqRecord,MySeqRecord,int, bool, bool, int]=None # originalUnsplitRecord, coonectedRecord,isLeft, isTop, splitPointIndex
		self.visualModel:tuple[int,int,int, int]=None # xStart, yStart, xStop, yStop. this could be used as unique id too probably. Bad style because we use visual position data in the model
		MySeqRecord.uniqueId+=1

	def __str__(self):
		collectDict = {}            
		# collectDict the attributes from the parent SeqRecord class
		collectDict['id'] = str(self.id)
		collectDict['name'] = str(self.name)
		collectDict['description'] = self.description+"\n"
		collectDict['size'] = str(len(self.seq))+"\n"
		collectDict['seq'] = self.seq._data.decode('ASCII')+"\n"
		featuresString="\n"
		for f in self.features:
			featuresString+="type:"+f.type+"   "
			if f.qualifiers.get("label"):
				featuresString+="location:"+str(f.location)+" text:"+f.qualifiers.get("label")[0]
				featuresString+="\n"
		collectDict['features'] = featuresString
		annotationsString="\n"
		#type id qualifier
		for key, value in self.annotations.items():
			annotationsString+=(f"{key}= {value}")        
			annotationsString+=", "
		annotationsString+="\n"
		collectDict['annotations'] = annotationsString        
		collectDict['dbxrefs'] =  str(self.dbxrefs)
		# collectDict additional attributes from MySeqRecord
		collectDict['\n*******'] = "*************************************************************************************\n"
		collectDict['isPrimer'] =  str(self.isPrimer)+"\n"
		collectDict['fiveTo3'] =  str(self.fiveTo3)+"\n"
		collectDict['xStartOffsetAsLetters'] =  str(self.xStartOffsetAsLetters)+"\n"
		if self.hybridizedToStrand:
			collectDict['hybridizedToStrand'] =  self.hybridizedToStrand.description+"-"+str("5 to 3 " if self.fiveTo3 else "3 to 5")
		if self.hybridizedToPrimer:
			collectDict['hybridizedToPrimer'] =  self.hybridizedToPrimer.description
		collectDict['uniqueId'] =  str(self.uniqueId)+"\n"
		collectDict['notAnnealedLocation'] =  str(self.notAnnealedLocation)+"\n"
		collectDict['singleStranded'] = str(self.singleStranded)+"\n"
		s=""
		for key in collectDict:
			s+=key+":"+collectDict[key]+" "
		#
		if self.loopInfo:
			s+=f"Loop info: Unsplit sequence:{self.loopInfo[0].seq},Other connected sequence:{self.loopInfo[1].seq}, left loop: {self.loopInfo[2]}, top loop:   {self.loopInfo[3]}\n"
		if self.visualModel:
			s+=f"Visual Model: Start y:{self.visualModel[0]}, end x:   {self.visualModel[2]}\n"
		return s

	def setNotAnnealedLocation(self, startStop:tuple):
		self.notAnnealedLocation: tuple[int, int]=startStop
		

	def getNotAnnealedLocation(self)->tuple:
		return self.notAnnealedLocation
		
	def __repr__(self):
		# Customize the representation of the item
		return "id"+str(self.id)+" description:"+self.description+" SingleStranded:"+str(self.singleStranded)+" primer:"+str(self.isPrimer)+" 5To3:"+str(self.fiveTo3) +" primer:"+str(self.isPrimer)+" Sequence:"+str( self.seq)[:70]


	def addFeature(self, start, end, strand:int, type, id, label:str):
		newFeature:SeqFeature=SeqFeature(SimpleLocation(start, end, strand=None), type=type, id=id, qualifiers={"label": [label]  })
		self.features.append(newFeature)   


	def removeFullSpanningFeatures(self): # that is usually where type=="source"
		for i in range(len(self.features) - 1, -1, -1):
			if self.toIgnore(self.features[i]):
					self.features.pop(i)	


	def removeFeaturesNotOnThisStrand(self): # that is usually where type=="source"
		for i in range(len(self.features) - 1, -1, -1):
			if self.fiveTo3 and self.features[i].strand != 1:
				self.features.pop(i)		
			if not self.fiveTo3 and self.features[i].strand != -1:
				self.features.pop(i)		

			
	def toIgnore(self,feature:SeqFeature)->bool: # usually for record.type=="source" that takes the full length of the sequence):
		if feature.location.start==0 and feature.location.end==len(self):
			return True
		return False

	def splitRecord(self, splitPointIndex)->MyT:
		# splits the record at splitPointIndex in left part and right part. 
		# calculates which one stays up and which one goes lower and becomes a primer
		# returns the top one first
		if splitPointIndex>=len(self.seq):
			return None, None
		seqString=MySeqRecord.seqToString(self.seq)
	
		seqR:Seq=  Seq(seqString[splitPointIndex:])
		seqL:Seq=  Seq(seqString[0:splitPointIndex])
		newFeaturesBefore, newFeaturesAfter=self.shiftFeaturesLocs(splitPointIndex)
		newRightSideRec=MySeqRecord(SeqRecord(seqR,id=self.id+"_t", name=self.name, description="truncated"+ self.description), singleStranded=True,fiveTo3=True,primer=False)
		newRightSideRec.features=newFeaturesAfter
		newRightSideRec.xStartOffsetAsLetters=self.xStartOffsetAsLetters+splitPointIndex
		newLeftSideRec=MySeqRecord(SeqRecord(seqL,id=self.id+"_l", name=self.name, description="looped"+ self.description), singleStranded=True,fiveTo3=True,primer=False)
		newLeftSideRec.features=newFeaturesBefore
		newLeftSideRec.xStartOffsetAsLetters=0#newUpper.xStartOffsetAsLetters+extraIndentForSecond is good too
		if  self.isLeft(splitPointIndex): # I split on the left side of visually shown sequence
			# leftSideDown=True
			newLeftSideRec.isPrimer=True
			return newRightSideRec, newLeftSideRec
		else:
			newRightSideRec.isPrimer=True
			return newLeftSideRec, newRightSideRec


	def isLeft(self, splitPointIndex):
		return splitPointIndex*2<len(self.seq)
	
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
		newFeaturesAfter=list()
		newFeaturesBefore=list()
		for f in self.features:
			f:SeqFeature			
			start=f.location.start
			end=f.location.end
			if start>byHowMuch:		# right side of the chopeed record that becomes the upper fragment		
				if end >byHowMuch:
					l=SimpleLocation(f.location.start-byHowMuch,f.location.end-byHowMuch)
				else:
					l=SimpleLocation(byHowMuch,f.location.end-byHowMuch)
				newFeature:SeqFeature=SeqFeature(l, type=f.type, id=f.id, qualifiers={"label": ["R "+f.qualifiers.get("label")[0]]  })
				newFeaturesAfter.append(newFeature)
			else:# left side of the chopped record that becomes the lower fragment	
				if end <=byHowMuch:
					l=SimpleLocation(start,end)
				else:
					l=SimpleLocation(start,byHowMuch)
				newFeature:SeqFeature=SeqFeature(l, type=f.type, id=f.id, qualifiers={"label": ["L "+f.qualifiers.get("label")[0]]  })
				newFeaturesBefore.append(newFeature)					
		return newFeaturesBefore, newFeaturesAfter
	
	def shiftBackwardsFeaturesLocs(self, lowerLoopLen)->list[SeqFeature]:
		newFeatures=list()
		for f in reversed(self.features):
			f:SeqFeature	

			end=lowerLoopLen-f.location.start
			start=lowerLoopLen-f.location.end
			if start>0:				
				if start >0:# fits the hole feature 
					l=SimpleLocation(start,end, strand=-1)
				# else:
				# 	l=SimpleLocation(byHowMuch,f.location.end-byHowMuch, strand=-1)
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
					l=SimpleLocation(f.location.start-byHowMuch,f.location.end-byHowMuch, strand=+1)
				else:
					l=SimpleLocation(byHowMuch,f.location.end-byHowMuch, strand=+1)
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
