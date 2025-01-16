
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition, UnknownPosition, Location

sequenceRecordList=[]
shrinkedSequenceRecordList=[]
loadedFileName:str

def appendSequenceRecord(newSequenceRecord:SeqRecord):
    print("Appended new record", newSequenceRecord)
    sequenceRecordList.append(newSequenceRecord)

def dumpModel(message):
    print(message,"Current Model", sequenceRecordList)
    for record in sequenceRecordList:
        print("Record id", record.id)
        print("Record name:", record.name)
        print("Record seq", record.seq)
        feature:SeqFeature
        for feature in record.features:
            print("Feature:",feature.id,  feature.location,feature.location.strand, feature.type )
            for q in feature.qualifiers:
                print("Qualifier:",q)
            loc:Location=feature.location
            print("Location:", loc.start, loc.end, loc.strand)


# 	for record in sequenceRecordIterator:
# 		#print("EMBL record name:%s length:%i Sequence: %s" % (record.id, len(record), record.seq))
# 	return 	sequenceRecordList[0].seq

# def appendSequence(newSequence:Seq):
# 	for record in sequenceRecordIterator:
# 		#print("EMBL record name:%s length:%i Sequence: %s" % (record.id, len(record), record.seq))
# 		sequenceRecordList.append(record)
# 	return 	sequenceRecordList[0].seq
# mod:Seq

