# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *              Fernando A. Teixeira (fernando.meireles@epfl.ch)
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
import csv
import re
from pwem.objects import AtomStruct, Sequence, SetOfSequences
from pyworkflow.constants import BETA
from pyworkflow.protocol import (params, 
                                LEVEL_ADVANCED,
                                GPU_LIST,
                                USE_GPU)
from pyworkflow.utils import Message
from pyworkflow.utils.path import makePath
from pyworkflow.object import Integer, String
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import AtomicStructHandler, fromCIFToPDB
from ..constants import CLUSTALO
from pwem.convert.sequence import alignClustalSequences
from Bio import SeqIO, AlignIO
from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1
from collections import Counter

from carbonara import Plugin


class AminoListParam(params.StringParam):
    """ list of aminoacids parameter
    """
    def __init__(self, **args):
        params.StringParam.__init__(self, **args)
        self.addValidator(AminoListValidator())


class AminoListValidator(params.Conditional):
    """ Validator for class AminoListParam(params.StringParam):
   """
    AMINO_LIST = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                  "M", "F", "P", "S", "T", "W", "Y", "V"]

    def __init__(self, error='Incorrect format for amino acid list param'):
        params.Conditional.__init__(self, error)

    def _condition(self, value):
        valid_pattern = re.compile(rf"^[{''.join(self.AMINO_LIST)},\s]+$")

        def validate_amino_sequence(seq):
            """Validate that the sequence contains only
               valid amino acids, commas, or spaces."""
            if not valid_pattern.match(seq):
                raise ValueError(f"Invalid character found in sequence: {seq}")
            # Extract amino acids only
            # this line is only for testing
            amino_acids = re.findall(rf"[{''.join(self.AMINO_LIST)}]", seq)
            print(f"Amino acids extracted: {amino_acids}")
            return True
        try:
            result = validate_amino_sequence(value)
            print(f"Validated amino acid sequence: {result}")
        except ValueError as e:
            print(f"Error: {e}")
            return False
        return True


class CarbonaraSamplingSequence(EMProtocol):
    """
    CARBonAra is an automatic method that generates multiple sequences
    able to fold in a certain structural configuration (atom structure used as input). 
    """
    _label = 'sampling sequence'
    _program = ""
    _version = ""

# -------------------------- DEFINE param functions ----------------------
    METHOD_OPTIONS = ['max', 'sampled']


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
        '''
        TODO: Restore these parameters when they work properly in the carbonara method
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
                     
        form.addParam('selectUnknounResidues', params.BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      label="Select specific RESIDUES for sampling?",
                      default=False,
                      help="Select 'Yes' if you want to select for sampling any specific residues "
                            "(one or more) of the atom structure.\n")
        
        form.addParam('selectUnknownStructureResidues', params.StringParam, 
                      expertLevel=LEVEL_ADVANCED, 
                      condition=('selectUnknounResidues==True'),
                      default=None, important=True,
                      label='Included residues',
                      help='"Use the wizard on the right to select the list of resides that will\n'
                           ' be overwritten to the two previous params and re-designed. New\n'
                           ' sequences will be only predicted for these selected positions and\n'
                           ' the rest of the sequence will be used as prior information for the\n'
                           ' prediction."')  
                    # Associate wizard to this parm  
        '''             
                     
        form.addParam('ignoreAminoacids', params.BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      label="Ignore specific AMINOACIDS from sampling?",
                      default=False,
                      help="Select 'Yes' if you want to exclude from sampling any specific aminoacids "
                            "(one or more) of the atom structure.\n")
        
        form.addParam('selectIgnoredAminoacids', AminoListParam,
                      expertLevel=LEVEL_ADVANCED,
                      condition=('ignoreAminoacids==True'),
                      default=None, important=True,
                      label='Ignored aminoacids',
                      help='"Write the list of aminoacids that will\n'
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
            args.extend(["--known_chains", chains.replace(" ", "")])
        '''
        TODO: Restore this parameter when it works in carbonara method

        # excluded residues        
        if self.selectKnoumResidues==True and self.selectKnownStructureResidues.get() is not None:
            knownResidues = str(self.selectKnownResidues.get())
            args.extend(["--known_positions", knownResidues])  
        
        # included residues  
        if self.selectUnknounResidues==True and self.selectUnknownStructureResidues.get() is not None:
            unknownResidues = str(self.selectUnknownStructureResidues.get())
            args.extend(["--unknown_positions", unknownResidues]) 
        '''
            
        # ignored aminoacid
        if self.ignoreAminoacids==True and self.selectIgnoredAminoacids.get() is not None:
            ignoredResidues = str(self.selectIgnoredAminoacids.get())
            args.extend(["--ignored_amino_acids", ignoredResidues.replace(" ", "")]) 
        
        # hetatm included
        if self.ignoreHetatm:
            args.extend(["--ignore_hetatm"]) 

        # water included
        if self.ignoreWater:
            args.extend(["--ignore_water"])

        # use GPU
        if not self.useGpu:
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
        output_scores_csv_file = os.path.join(self.dir_path, "sorted_scores.csv")
        self.extract_scores_from_fasta(self.subdir_path, output_scores_csv_file)

    def createOutputStep(self):
        """Register sequences generated"""

        setSeq = SetOfSequences()
        outputSequences = setSeq.create(outputPath=self._getPath())

        sequences = []
        # all the sampled sequences without duplicates
        filepath_list = self.sel_files_seq_unique(self.subdir_path)

        # check if .fasta files exist before registering
        # in a folder call sequences
        # create output: sequences
        for filename in sorted(os.listdir(self.subdir_path)):
            if filename.endswith(".fasta"):
                path = os.path.join(self.subdir_path, filename)
                sequence = self.extract_sequence_from_file(path)
                label= os.path.splitext(filename)[0]

                if path in filepath_list:
                    seq = Sequence()
                    seq.setSeqName(path)
                    seq.setId(label)
                    seq.setSequence(sequence)
                    seq.setObjLabel(label)
                    keyword = filename.split(".fasta")[0]
                    kwargs = {keyword: seq}
                    self._defineOutputs(**kwargs)

                    # create output: setOfSequences
                    sequences.append(seq)
        
        for seq in sequences:
            outputSequences.append(seq)
        
        self._defineOutputs(outputSequences=outputSequences)

        # total number of unique sequences
        num_sequences = len(sequences)
        print(f"Total number of unique sequences: {num_sequences}")
   

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
            summary.append("%s *sequences* predicted " % counter)
            summary.append("")
            filepath_list = self.sel_files_seq_unique(subdir_path)
            number_unique_sequences = len (filepath_list)
            summary.append("%s *unique sequences* predicted " % number_unique_sequences)
            summary.append("")
            summary.append("Alignment generated with Clustal Omega and saved in:")
            
            clustal_aln_path = os.path.join(dir_path, "clustal.aln")
            summary.append(" %s " % clustal_aln_path)
            summary.append("")

            # Showing ClustalOmega alignment symmary with a monospaced font 

            clustal_aln_summmary_path = os.path.join(dir_path, "clustal_summary.aln")
            try:
                with open(clustal_aln_summmary_path, "r") as f:
                    for line in f:
                        if line == "\n":
                            continue
                        summary.append("'''" + line.rstrip('\n') + "'''\n")

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

            #include all the sampled sequences next without duplicates
            filepath_list = self.sel_files_seq_unique(folder_path)
            
            for filename in os.listdir(folder_path):
                if filename.endswith('.fasta'):
                    filepath = os.path.join(folder_path, filename)
                    if filepath in filepath_list:
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
    
    def extract_scores_from_fasta(self, folder_path, output_scores_csv_file):
        """
        Generate a CSV file with columns: file_name, sequence, score.
        Avoids duplicate sequences and sorts entries by score (highest first).
        """
        rows = []
        
        #include all the sampled sequences without duplicates
        filepath_list = self.sel_files_seq_unique(folder_path)

        for filename in sorted(os.listdir(folder_path)):
            if filename.endswith(".fasta"):
                filepath = os.path.join(folder_path, filename)

                if filepath in filepath_list:
                    for record in SeqIO.parse(filepath, "fasta"):
                        seq_str = str(record.seq)

                        #extract score from header
                        header = record.description
                        score = None
                        for part in header.split(","):
                            if "score" in part:
                                score = float(part.split("=")[1].strip())
                                rows.append([filename, seq_str, score])
                                break

        # sort rows by score descending
        rows.sort(key=lambda x: x[2], reverse=True)
        
        # write to CSV
        with open(output_scores_csv_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["file_name", "sequence", "CARBonAra score"]) # header row
            writer.writerows(rows)

        print(f"CSV file generated: {output_scores_csv_file} with {len(rows)} unique sequences")

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

    def extract_sequence_from_file(self, filepath):
        sequence_lines = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('>'):
                    continue  
                sequence_lines.append(line)
        return ''.join(sequence_lines)
    
    def sel_files_seq_unique(self, folder_path):
        # method to select CARBonAra output sequence files avoiding sequence duplicates
        seq_unique = []
        filepaths = []
        for filename in sorted(os.listdir(folder_path)):
            if filename.endswith(".fasta"):
                filepath = os.path.join(folder_path, filename)

                for record in SeqIO.parse(filepath, "fasta"):
                    seq_str = str(record.seq)
                    if seq_str not in seq_unique:
                        seq_unique.append(seq_str)
                        filepaths.append(filepath)

        return filepaths


