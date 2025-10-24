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

import os
import glob
from pwem.objects import AtomStruct, Sequence, SetOfSequences
from pyworkflow.constants import BETA
from pyworkflow.protocol import (params, 
                                LEVEL_ADVANCED,
                                GPU_LIST,
                                USE_GPU)
from pyworkflow.utils import Message
from pyworkflow.object import Integer
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import AtomicStructHandler, fromCIFToPDB
from ..constants import CLUSTALO
from pwem.convert.sequence import alignClustalSequences
from Bio import SeqIO, AlignIO
from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1
from collections import Counter

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
    METHOD_OPTIONS=['max', 'sampled']
    AMINO_LIST=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

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
                      default=0.5,
                      label='Ratio of prior information', important=True,
                      help='Ratio of sequence imprint for sampling'
                           'Value between 0.0 and 1.0. Ratio of the prediction to be used as\n'
                           ' prior information. To get higher variability in the sequence you\n'
                           ' have to choose values close to 1.0. All the positions will contain\n'
                           ' prior information to bias the prediction')
    
        form.addParam('bSampled', params.EnumParam, 
                      choices=self.METHOD_OPTIONS,
                      display=params.EnumParam.DISPLAY_LIST,
                      default=1,
                      expertLevel=LEVEL_ADVANCED, allowsNull=False,
                      label='Sampling method',
                      help='"sampled" means sampling from probability. It will result in higher\n'
                           ' diversity while still maintaining reasonable sequence confidence.\n'
                           ' "max" means sampling with maximum confidence. It will result in\n'
                           ' high-confidence predictions but low diversity.')
        
        form.addParam('selectChains', params.BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      label="Exclude any CHAINS from sampling?",
                      default=False,
                      help="Select 'Yes' if you want to exclude from sampling any specific chain "
                            "(one or more) of the atom structure.\n")
        
        form.addParam('selectStructureChains', params.StringParam,
                      expertLevel=LEVEL_ADVANCED, 
                      condition=('selectChains==True'),
                      default=None, important=True,
                      label='Excluded chains',
                      help='"Use the wizard on the right to select the list of known chains for\n'
                           ' which no new sequences will be predicted. Moreover, the sequence\n'
                           ' information of the selected chains will be used as prior information\n'
                           ' for the prediction. "')
                    # Associate wizard to this param

        form.addParam('selectKnoumResidues', params.BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      label="Exclude any RESIDUES from sampling?",
                      default=False,
                      help="Select 'Yes' if you want to exclude from sampling any specific residues "
                            "(one or more) of the atom structure.\n")
                      
        form.addParam('selectKnownStructureResidues', params.StringParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      condition=('selectKnoumResidues==True'),
                      default=None, important=True,
                      label='Excluded residues',
                      help='"Use the wizard on the right to select the list of known resides for\n'
                           ' which no new sequences will be predicted. Moreover, the sequence\n'
                           ' information of the selected positions will be used as prior\n'
                           ' information for the prediction. "')  
        # Associate wizard to this parm   
                      
        form.addParam('selectUnknoumResidues', params.BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      label="Select specific RESIDUES for sampling?",
                      default=False,
                      help="Select 'Yes' if you want to select for sampling any specific residues "
                            "(one or more) of the atom structure.\n")
        
        form.addParam('selectUnknownStructureResidues', params.StringParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      condition=('selectUnknoumResidues==True'),
                      default=None, important=True,
                      label='Included residues',
                      help='"Use the wizard on the right to select the list of resides that will\n'
                           ' be overwritten to the two previous params and re-designed. New\n'
                           ' sequences will be only predicted for these selected positions and\n'
                           ' the rest of the sequence will be used as prior information for the\n'
                           ' prediction."')  
                    # Associate wizard to this parm  
                      
        form.addParam('ignoreAminoacids', params.BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore specific AMINOACIDS from sampling?",
                      default=False,
                      help="Select 'Yes' if you want to exclude from sampling any specific aminoacids "
                            "(one or more) of the atom structure.\n")
        
        form.addParam('selectIgnoredAminoacids', params.StringParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      choices=self.AMINO_LIST, 
                      condition=('ignoreAminoacids==True'),
                      default=None, important=True,
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

        # get atom structure PDBx/mmCIF file name

        fileName = self.atomStruct.get().getFileName()
        self.atomStructName = os.path.abspath(fileName)
        
        # convert CIF to PDB file
        if self.atomStructName.endswith('.cif'):
            atomStructReName = os.path.abspath(
                self._getExtraPath(os.path.basename(
                    os.path.splitext(self.atomStructName)[0] + ".pdb")))
            fromCIFToPDB(self.atomStructName, atomStructReName, log = self._log)
            self.atomStructName =  atomStructReName
        
        # get a dict of chains with residues
        h = AtomicStructHandler()
        h.read(self.atomStructName) 
        h.getStructure()
        self.listOfChains, self.listOfResidues = h.getModelsChains()
        print("self.listOfChains: ", self.listOfChains)

    def processStep(self):
    
        args = []
        
        # number of atom structures
        args.extend(["--num_sequences", str(self.numSamples.get())])
        
        # ratio of prior information
        args.extend(["--imprint_ratio", str(self.imprintRadio.get())])
                    
        # sampling method  
        args.extend(["--sampling_method", self.METHOD_OPTIONS[self.bSampled.get()]])

        # excluded chains        
        if self.selectChains==True and self.selectStructureChains.get() is not None:
            chains = str(self.selectStructureChains.get())
            args.extend(["--known_chains", chains])

        # excluded residues
        if self.selectKnoumResidues==True and self.selectKnownStructureResidues.get() is not None:
            knownResidues = str(self.selectKnownResidues.get())
            args.extend(["--known_positions", knownResidues])  

        # included residues  
        if self.selectUnknoumResidues==True and self.selectUnknownStructureResidues.get() is not None:
            unknownResidues = str(self.selectUnknownStructureResidues.get())
            args.extend(["--unknown_positions", unknownResidues]) 

        # ignored aminoacid
        if self.ignoreAminoacids==True and self.selectIgnoredAminoacids.get() is not None:
            ignoredResidues = str(self.selectIgnoredAminoacids.get())
            args.extend(["--ignored_amino_acids", ignoredResidues]) 
        
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
            os.environ["CUDA_VISIBLE_DEVICES"] = ("%s" % self.getGpuList()[0])  # selected GPU
            print("I'm using GPU " + ("%s" % self.getGpuList()[0]))

        # pdb_filepath
        args.extend([self.atomStructName])

        # output_dir, folder called carbonara_results
        self.dir_path = os.path.abspath(self._getExtraPath())
        subdirectory = "carbonara_results"
        self.subdir_path = os.path.join(self.dir_path, subdirectory)
        args.extend([self.subdir_path])

        args_str = ' '.join(args)

        # Call carbonara:
        self.runJob(Plugin.getCarbonaraCmd(), args_str)

        # Align output sequences and save the alignment
        self.clustalOAlignment()

        # Summary of alignment with consensus
        summ_align_file = os.path.abspath(self._getExtraPath("clustal_summary.aln")) 
        self.filter_file_lines(self.outFile, summ_align_file)

        # Summarize sequence scores in a file
        output_scores_file = os.path.join(self.dir_path, "sorted_scores.txt")
        self.extract_scores_from_fasta(self.subdir_path, output_scores_file)


    def createOutputStep(self):
        """Register sequences generated"""
        # check if .fasta files exist before registering
        # in a folder call sequences

        # TODO: save data as a set of sequences?
        for filename in sorted(os.listdir(self.subdir_path)):
            if filename.endswith(".fasta"):
                path = os.path.join(self.subdir_path, filename)
                seq = Sequence()
                seq.setName(path)
                keyword = filename.split(".fasta")[0]
                kwargs = {keyword: seq}
                self._defineOutputs(**kwargs)

        

    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        errors = []

        gpus = self.getGpuList()

        if len(gpus) > 1:
            errors.append('Only one GPU can be used.')

        # Check that CLUSTALO program exists
        if not (self.is_tool(CLUSTALO)):
            errors.append("Clustal-omega program missing.\n "
                          "You need it to run this program.\n"
                          "Please install Clustal-omega:\n"
                          "     sudo apt-get install clustalo\n")

        return errors

    def _summary(self):

        summary = []
        if self.getOutputsSize() > 0:
            dir_path = os.path.abspath(self._getExtraPath())
            subdirectory = "carbonara_results"
            subdir_path = os.path.join(dir_path, subdirectory)
            counter = 0
            for filename in sorted(os.listdir(subdir_path)):
                if filename.endswith(".fasta"):
                    counter += 1       
            summary.append("%s sequences predicted" % counter)
            summary.append("" )
            summary.append("Alignment generated with Clustal Omega and saved in:" )
            clustal_aln_path = os.path.join(dir_path, "clustal.aln")
            summary.append(" %s " % clustal_aln_path)
            summary.append("" )

            # This part only makes sense if a monospaced font like courier, consolas, monaco, lucida console
            # otherwise I have to modify it to show only the consensus line
            clustal_aln_summmary_path = os.path.join(dir_path, "clustal_summary.aln")
            try:
                with open(clustal_aln_summmary_path, "r") as f:
                    content = f.read()
                summary.append(content )
            except FileNotFoundError:
                summary.append(f"File not found: {clustal_aln_summmary_path}")
            except Exception as e:
                summary.append(f"Error reading file: {e}")
            
        else:
            summary.append(Message.TEXT_NO_OUTPUT_FILES)
        return summary
    
    def _citations(self):
        return ['Krapp2024-dw']
    
    def clustalOAlignment(self):

        # generate file with fasta sequences
        inFile = self.merge_fasta_to_align(self.subdir_path) 

        # generate output file and run clustal
        self.outFile = os.path.abspath(self._getExtraPath("clustal.aln")) 
        cline = alignClustalSequences(inFile, self.outFile)  

        # new clustalo command line
        # clustalo -i input.fasta -o aligned.aln --outfmt=clustal --force 
        args = ' --outfmt=clustal --force '
        # cline is a biophythom object let us
        # convert it to a str before using it in scipion
        # TODO: Solve str bug in chimera modeller
        self.runJob(str(cline), args)
 
        
    def merge_fasta_to_align(self, folder_path):
        inFile =  os.path.abspath(self._getExtraPath("merged.fasta")) 
        self.wholeSequence = self.extract_protein_sequences(self.atomStructName)
        with open(inFile, 'w') as out_f:

            # include starting sequence as first sequence
            name = os.path.basename(
                    os.path.splitext(self.atomStructName)[0])
            out_f.write(f'>{name}\n{str(self.wholeSequence)}\n')

            #include all the sampled sequences next
            for filename in os.listdir(folder_path):
                if filename.endswith('.fasta'):
                    filepath = os.path.join(folder_path, filename)
                    name = os.path.splitext(filename)[0]
                    for record in SeqIO.parse(filepath, 'fasta'):
                        out_f.write(f'>{name}\n{str(record.seq)}\n')
        return inFile

    def is_tool(self, name):
        """Check whether `name` is on PATH."""
        from shutil import which
        return which(name) is not None
    
    def extract_protein_sequences(self, pdb_file):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file)

        chain_sequences = []

        for model in structure:
            for chain in model:
                residues = [res for res in chain if is_aa(res, standard=True)]
                if residues:
                    sequence = ''.join(seq1(res.get_resname()) for res in residues)
                    chain_sequences.append(sequence)
            break  # Only process the first model

        return ':'.join(chain_sequences)
    
    def extract_scores_from_fasta(self, folder_path, output_scores_file):
        entries = []

        # Find all .fasta files in the folder
        for filepath in glob.glob(os.path.join(folder_path, "*.fasta")):
            filename = os.path.basename(filepath)
            name_without_ext = os.path.splitext(filename)[0]

            with open(filepath, "r") as f:
                header = f.readline().strip()
                if header.startswith(">"):
                    # Extract score from header line
                    parts = header[1:].split(",")
                    score_part = next((p for p in parts if "score=" in p), None)
                    if score_part:
                        try:
                            score = float(score_part.split("=")[1])
                            entries.append((name_without_ext, score))
                        except ValueError:
                            print(f"Invalid score in {filename}")

        # Sort entries by score descending
        entries.sort(key=lambda x: x[1], reverse=True)

        # Write to output file
        with open(output_scores_file, "w") as out:
            for name, score in entries:
                out.write(f"{name}\t{score:.5f}\n")

    def filter_file_lines(self, filepath, newfilepath):
        first_label = None
        filtered_lines = []

        try:
            with open(filepath, "r") as f:
                for line in f:
                    stripped = line.strip()
                    # Keep header and empty lines
                    if not stripped or stripped.startswith("CLUSTAL"):
                        filtered_lines.append(line)
                        continue 

                    # Keep consensus lines
                    if stripped[0] in "*:." or line.startswith(" "):
                        filtered_lines.append(line)
                        continue

                    # Process sequence lines
                    parts = line.rstrip().split()
                    if len(parts) == 2:
                        label, seq = parts
                        if first_label is None:
                            first_label = label
                            filtered_lines.append(line)
                        elif first_label and label == first_label:
                            filtered_lines.append(line)
                    # Skip other labels

            # New file with valid lines
            with open(newfilepath, "w") as f:
                f.writelines(filtered_lines)

        except Exception as e:
            print(f"Error processing file: {e}")
    

