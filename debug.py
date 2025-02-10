import tkinter as tk
from copy import deepcopy
import math
from multiprocessing import log_to_stderr
from tkinter import Canvas, Button
import tkinter as tk
from tkinter import filedialog, messagebox, Menu
#
import sys
from pathlib import Path
from turtle import left
from typing import no_type_check  
from PIL import ImageGrab
#
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import time
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition, UnknownPosition, Location
# mine
import gl
from util import loadSequencesFile
from gl import *
from preferences import Preferences

