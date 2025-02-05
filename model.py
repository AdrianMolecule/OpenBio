from myseqrecord import MySeqRecord


class Model:
    # Class attribute of type Model
    modelInstance: "Model" = None  # Initially set to None

    def __init__(self, loadedFileName:str,sequenceRecordList: list = []):
        # Instance attributes
        self.sequenceRecordList: list[MySeqRecord] = sequenceRecordList
        self.loadedFileName: str = loadedFileName
        
    def dumpModel(self, optionalMessage=None) -> dict:
        """Method to display the current state of the model."""
        m= {"where":optionalMessage,
            "loadedFileName": self.loadedFileName,
            "sequenceRecordList": self.sequenceRecordList,
        }
        print( m)
        
    @classmethod
    def setModel(cls, instance: "Model") -> None:
        """Class method to set the class attribute modelInstance."""
        if isinstance(instance, Model): 
            cls.modelInstance = instance
        else:
            raise ValueError("The instance must be of type Model.")

# Example of using the Model class
if __name__ == "__main__":
    # Create an instance of Model
    loadedFileName = "data.txt"
    # Set some instance attributes
    loadedFileName = "data.txt"
    sequenceRecordList = ["record1", "record2", "record3"]
    shrinkedSequenceRecordList = ["record1", "record2"]
    myModel = Model(loadedFileName,sequenceRecordList)
    # Set the class attribute to the instance
    Model.setModel(myModel)
    # Call dumpModel to get the current state
    modelDump = myModel.modelInstance.dumpModel()
    print(Model.modelInstance)
    print(modelDump)
