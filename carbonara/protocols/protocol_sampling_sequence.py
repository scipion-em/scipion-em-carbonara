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


"""
This module allows to run CARBonAra method to generate multiples sequences based on a certain atom structure"
"""
from enum import Enum

from pwem.objects import AtomStruct, Sequence, SetOfSequences
from pyworkflow.constants import BETA
from pyworkflow.protocol import params, LEVEL_ADVANCED
from pyworkflow.utils import Message
from pyworkflow.object import Integer
from pwem.protocols import EMProtocol
from pyworkflow.protocol import GPU_LIST, USE_GPU

from carbonara import Plugin


class outputs(Enum):
    count = Integer


class ProtCarbonara(EMProtocol):
    """
    CARBonAra is an automatic method that generates sequences using an atom structure as input. 
    """
    _label = 'carbonara'
    _devStatus = BETA
    _possibleOutputs = {
    	OUTPUT_NAME: AtomStruct
    }

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputAtomStructure', params.PointerParam,
                      label='Atom Structure', important=True,
                      help='Add here the atom structure for which you want to generate\n' 
                           ' sequences able to fold in agreement with the architecture\n' 
                           ' of that scaffold structure. The program automatically removes\n'
                           ' hydrogens and side chains from the scaffold. The default is to\n' 
                           ' keep HETATM and water, but this behavior can be toggled.') 

        form.addParam('numSample', params.IntParam,
                      validators=[params.Positive],
                      default=100,
                      label='Number of Sequences', important=True,
                      help='Number of sequences to sample.')
                      
        form.addParam('imprintRadio', params.FloatParam,
                      validators=[params.Positive],
                      default=0.5,
                      label='Ratio of prior information', important=True,
                      help='Value between 0.0 and 1.0. Ratio of the prediction to be used as\n'
                           ' prior information. To get higher variability in the sequence you\n'
                           ' have to choose values close to 1.0. All the positions will contain\n'
                           ' prior information to bias the prediction')

        form.addParam('BSampled', params.BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label='Maximum diversity?',
                      help='"Yes" means sampling from probability. It will result in higher\n'
                           ' diversity while still maintaining reasonable sequence confidence.\n'
                           ' "No" means sampling with maximum confidence. It will result in\n'
                           ' high-confidence predictions but low diversity.')
        
        form.addParam('selectStructureChain', params.MultiPointerParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      pointerClass="Chain", allowsNull=True,
                      label='Excluded chains',
                      help="Use the wizard on the right to select the list of known chains for\n'
                           ' which no new sequences will be predicted. Moreover, the sequence\n'
                           ' information of the selected chains will be used as prior information\n'
                           ' for the prediction. "')   
                      
        form.addParam('selectKnownResidues', params.MultiPointerParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      pointerClass="Residue", allowsNull=True,
                      label='Excluded residues',
                      help="Use the wizard on the right to select the list of known resides for\n'
                           ' which no new sequences will be predicted. Moreover, the sequence\n
                           ' information of the selected positions will be used as prior\n'
                           ' information for the prediction. "')     
                      
        form.addParam('selectUnknownResidues', params.MultiPointerParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      pointerClass="Residue", allowsNull=True,
                      label='Included residues',
                      help="Use the wizard on the right to select the list of resides that will\n'
                           ' be overwritten to the two previous params and re-designed. New\n'
                           ' sequences will be only predicted for these selected positions and\n'
                           ' the rest of the sequence will be used as prior information for the\n'
                           ' prediction."')    
                      
        form.addParam('selectIgnoredResidues', params.MultiPointerParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      pointerClass="Residue", allowsNull=True,
                      label='Ignored residues',
                      help="Use the wizard on the right to select the list of resides that will\n'
                           ' be completely ignored for the sequence sampling. The prior\n'
                           ' information and sequences generated will not contain the selected\n'
                           ' amino acids."') 
                      
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
                       label="Use GPU for execution?",
                       help="This protocol has both CPU and GPU implementation. "
                            "Select the one you want to use.")

        form.addHidden(GPU_LIST, params.StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU ID (single one)",
                       help="GPU device to be used")

    # --------------------------- INSERT steps functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.predictStep, needsGPU=True)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions ------------------------------
    def greetingsStep(self):
        # say what the parameter says!!
        for time in range(0, self.times.get()):
            print(self.message)

    def createOutputStep(self):
        # register how many times the message has been printed
        # Now count will be an accumulated value
        timesPrinted = Integer(self.times.get() + self.previousCount.get())

        self._defineOutputs(**{outputs.count.name: timesPrinted})
        self._defineSourceRelation(self.message, timesPrinted)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []

        if self.times > 20:
            errors.append("Cannot do more than 20 times.")

        return errors

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []

        if self.isFinished():
            summary.append(f"This protocol has printed *{self.message}* {self.times} times.")
        return summary

    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append(f"{self.message} has been printed in this run {self.times} times.")
            if self.previousCount.hasPointer():
                methods.append("Accumulated count from previous runs were %i."
                               " In total, %s messages has been printed."
                               % (self.previousCount, self.count))
        return methods
