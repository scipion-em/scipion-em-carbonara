# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * National Center for Biotechnology (CNB-CSIC)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from enum import Enum
import os 
import subprocess
from pwem.objects import AtomStruct, Sequence, SetOfSequences
from pwem.convert.atom_struct import fromCIFToPDB
from pyworkflow.constants import BETA
from pyworkflow.protocol import (params, 
                                LEVEL_ADVANCED,
                                GPU_LIST,
                                USE_GPU)
from pyworkflow.utils import Message
from pyworkflow.object import Integer
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import AtomicStructHandler


from carbonara import Plugin

class CarbonaraSamplingSequence(EMProtocol):
    """
    CARBonAra is an automatic method that generates multiple sequences
    able to fold in a certain structural configuration (atom structure used as input). 
    """
    _label = 'sampling sequence'
    _program = ""
    _version = ""

# -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Defining the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('atomStruct', params.PointerParam,
                      label='Atom Structure', important=True,
                      pointerClass="AtomStruct", allowsNull=False,
                      help='Add here the atom structure for which you want to generate\n' 
                           ' sequences able to fold in agreement with the architecture\n' 
                           ' of that scaffold structure. The program automatically removes\n'
                           ' hydrogens and side chains from the scaffold. The default is to\n' 
                           ' keep HETATM and water, but this behavior can be toggled.') 

        form.addParam('numSamples', params.IntParam,
                      validators=[params.Positive],
                      default=100,
                      label='Number of Sequences', important=True,
                      help='Number of sequences to generate.')
                      
        form.addParam('imprintRadio', params.FloatParam,
                      validators=[params.Positive],
                      default=0.5,
                      label='Ratio of prior information', important=True,
                      help='Ratio of sequence imprint for sampling'
                           'Value between 0.0 and 1.0. Ratio of the prediction to be used as\n'
                           ' prior information. To get higher variability in the sequence you\n'
                           ' have to choose values close to 1.0. All the positions will contain\n'
                           ' prior information to bias the prediction')

        form.addParam('bSampled', params.EnumParam, 
                      choices=['max', 'sampled'],default=1,
                      expertLevel=LEVEL_ADVANCED, allowsNull=False,
                      label='Sampling method',
                      help='"sampled" means sampling from probability. It will result in higher\n'
                           ' diversity while still maintaining reasonable sequence confidence.\n'
                           ' "max" means sampling with maximum confidence. It will result in\n'
                           ' high-confidence predictions but low diversity.')
        
        form.addParam('selectStructureChains', params.EnumParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      choices=[], allowsNull=True,
                      label='Excluded chains',
                      help='"Use the wizard on the right to select the list of known chains for\n'
                           ' which no new sequences will be predicted. Moreover, the sequence\n'
                           ' information of the selected chains will be used as prior information\n'
                           ' for the prediction. "')   
                    # Associate wizard to this parm
                      
        form.addParam('selectKnownResidues', params.EnumParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      choices=[], allowsNull=True,
                      label='Excluded residues',
                      help='"Use the wizard on the right to select the list of known resides for\n'
                           ' which no new sequences will be predicted. Moreover, the sequence\n'
                           ' information of the selected positions will be used as prior\n'
                           ' information for the prediction. "')  
                    # Associate wizard to this parm   
                      
        form.addParam('selectUnknownResidues', params.EnumParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      choices=[], allowsNull=True,
                      label='Included residues',
                      help='"Use the wizard on the right to select the list of resides that will\n'
                           ' be overwritten to the two previous params and re-designed. New\n'
                           ' sequences will be only predicted for these selected positions and\n'
                           ' the rest of the sequence will be used as prior information for the\n'
                           ' prediction."')  
                    # Associate wizard to this parm  
                      
        form.addParam('selectIgnoredAminoacids', params.EnumParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      choices=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"], 
                      allowsNull=True,
                      label='Ignored aminoacids',
                      help='"Use the wizard on the right to select the list of aminoacids that will\n'
                           ' be completely ignored for the sequence sampling. The prior\n'
                           ' information and sequences generated will not contain the selected\n'
                           ' amino acids."') 
                    # Associate wizard to this parm
                      
        form.addParam('ignoreHetatm', params.BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Ignore HETATM?',
                      help='By default every atoms included in the input file are used for the\n'
                           ' prediction, such as ligands, lipids, ions, and water.')    
                      
        form.addParam('ignoreWater', params.BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Ignore water molecules?',
                      help='By default water molecules included in the input file are used for\n'
                           ' the prediction.')         
                                
        form.addHidden(USE_GPU, params.BooleanParam, default=True,
                       expertLevel=LEVEL_ADVANCED,
                       label="Use GPU for execution?",
                       help="This protocol has both CPU and GPU implementation. "
                            "Select the one you want to use.")

        form.addHidden(GPU_LIST, params.StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       condition=('USE_GPU==True'),
                       label="Choose GPU ID (single one)",
                       help="GPU device to be used")

    # --------------------------- INSERT steps functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep(self.preRequisitesStep, needsGPU=False)
        self._insertFunctionStep(self.processStep, needsGPU=True)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions ------------------------------
    def preRequisitesStep(self):

        # get atom structure file name
        self.atomStructure = self.atomStruct.get()
        self.atomStructName = os.path.abspath(self.atomStructure.getFileName())
        self.inputStructure

        fromCIFToPDB(self.atomStructName, self.inputStructure)

        structureHandler = AtomicStructHandler()
        structureHandler.read(self.inputStructure) 
        structureHandler.getStructure()
        self.listOfChains, self.listOfResidues = structureHandler.getModelsChains()
        # returns a dict of chains with residues

    def processStep(self):
    
        args = []
        
        # number of atom structures
        args.extend(["--num_sequences", self.numSample])
        
        # ratio of prior information
        args.extend(["--imprint_ratio", self.imprintRadio])
                    
        # sampling method  
        args.extend(["--sampling_method", self.bSampled])

        # excluded chains        
        if len(self.selectStructureChains) > 0:
            args.extend(["--known_chains", self.selectStructureChains])

        # excluded residues
        if len(self.selectKnownResidues) > 0:
            args.extend(["--known_positions", self.selectKnownResidues])  

        # included residues  
        if len(self.selectUnknownResidues) > 0:
            args.extend(["--unknown_positions", self.selectUnknownResidues]) 

        # ignored aminoacid
        if len(self.selectIgnoredAminoacids) > 0:
            args.extend(["--ignored_amino_acids", self.selectIgnoredAminoacids]) 
        
        # hetatm included
        if self.ignoreHetatm:
            args.extend(["--ignore_hetatm"]) 

        # water included
        if self.ignoreWater:
            args.extend(["--ignore_water"])

        # use GPU
        if not USE_GPU:
            args.extend(["--device", "cpu"])
        else:
            args.extend(["--device", "cuda"])
            os.environ["CUDA_VISIBLE_DEVICES"] = self.getGpuList()[0]  # selected GPU
            print("I'm using GPU " + self.getGpuList()[0])

        # pdb_filepath
        args.extend([self.inputStructure])

        # output_dir
        args.extend([self._getExtraPath()])

        # Call carbonara:
        self.runJob(Plugin.getCarbonaraCmd(), args)

        
    def createOutputStep(self):
        """Register sequences generated"""
        # check if .fasta files exist before registering
        directory = self._getExtraPath()

        for filename in sorted(os.listdir(directory)):
            if filename.endswith(".fasta"):
                path = os.path.join(directory, filename)
                seq = Sequence()
                seq.setSeqNames(path)
                keyword = filename.split(".fasta")[0]
                kwargs = {keyword: seq}
                self._defineOutputs(**kwargs)

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        errors = []

        gpus = self.getGpuList()

        if len(gpus) > 1:
            errors.append('Only one GPU can be used.')

        return errors

    def _summary(self):
        # Think on how to update this summary with created PDB
        summary = []
        if self.getOutputsSize() > 0:
            directory = self._getExtraPath()
            counter = 0
            for filename in sorted(os.listdir(directory)):
                if filename.endswith(".fasta"):
                    counter += 1       
            summary.append("%s sequence predicted" % counter)
        else:
            summary.append(Message.TEXT_NO_OUTPUT_FILES)
        return summary
    
    def _citations(self):
        return ['Krapp2024-dw']