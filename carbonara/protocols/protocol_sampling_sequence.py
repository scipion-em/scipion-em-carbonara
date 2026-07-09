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

"""
Protocol implementation for CARBonAra sequence sampling.

This module defines CarbonaraSamplingSequence, the main Scipion protocol
that wraps the CARBonAra protein sequence design tool. The protocol:

  1. preRequisitesStep: converts CIF to PDB if needed, parses chain info
  2. processStep: builds CLI args, runs CARBonAra, performs ClustalO
     alignment, optionally runs ColabFold/AlphaFold scoring
  3. createOutputStep: registers generated sequences as Scipion outputs

Also defines AminoListParam and AminoListValidator for validating
comma-separated amino acid lists in the protocol form.
"""

import os
import csv
import re
import json
import subprocess
from pwem.objects import AtomStruct, Sequence, SetOfSequences
from pyworkflow.protocol import (params,
                                 LEVEL_ADVANCED,
                                 GPU_LIST,
                                 USE_GPU)
from pyworkflow.utils import Message, greenStr
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import AtomicStructHandler, fromCIFToPDB
from ..constants import CLUSTALO, conda_env
from pwem.convert.sequence import alignClustalSequences
from Bio import SeqIO
from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1

from carbonara import Plugin


class AminoListParam(params.StringParam):
    """String parameter that validates a comma/space-separated list of
    single-letter amino acid codes (e.g. 'A, R, N, D, C')."""
    def __init__(self, **args):
        params.StringParam.__init__(self, **args)
        self.addValidator(AminoListValidator())


class AminoListValidator(params.Conditional):
    """Validator for AminoListParam.

    Checks that the input string contains only valid single-letter amino
    acid codes, commas, and whitespace. Uses a regex built from the
    canonical 20 standard amino acids.
    """
    AMINO_LIST = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                  "M", "F", "P", "S", "T", "W", "Y", "V"]

    def __init__(self, error='Incorrect format for amino acid list param'):
        params.Conditional.__init__(self, error)

    def _condition(self, value):
        """Return True if value contains only valid amino acid letters,
        commas, and whitespace."""
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
        try:
            result = validate_amino_sequence(value)
            self._log.info(f"Validated amino acid sequence: {result}")
        except ValueError as e:
            self._log.error(f"Validation failed: {e}")
            return False
        return True


class CarbonaraSamplingSequence(EMProtocol):
    """
    CARBonAra is an automatic method that generates multiple sequences
    able to fold in a certain structural configuration (atom structure used as input).

    AI Generated:

    Carbonara Sampling Sequence (CarbonaraSamplingSequence) — User Manual
        Overview

        The Carbonara Sampling Sequence protocol is designed to generate new
        protein sequences that are compatible with a given structural scaffold.
        Starting from an experimentally determined or computationally derived
        atomic structure, the protocol predicts alternative amino acid sequences
        that are expected to preserve the overall fold and structural organization
        of the original protein architecture.

        This type of workflow is especially useful in protein engineering,
        binder optimization, sequence diversification, and exploratory design
        studies where the user wants to investigate which sequence variants may
        remain structurally compatible with a target fold. The protocol combines
        structural information with probabilistic sequence sampling strategies
        to generate multiple candidate sequences with varying levels of diversity
        and confidence.

        Inputs and Structural Context

        The main input is an atomic structure in PDB or mmCIF format. The protocol
        interprets this structure as the scaffold that defines the spatial and
        biochemical imprintRationment guiding sequence generation. During preprocessing,
        hydrogens and side-chain information are simplified to focus the prediction
        on the structural backbone and overall architecture.

        The protocol supports both single-chain proteins and multimeric complexes.
        In complexes, users may preserve selected chains as fixed structural
        context while redesigning only specific partners. This is particularly
        relevant in binder design workflows, interface optimization, or studies
        of protein-protein recognition.

        Sequence Sampling Strategy

        Sequence generation is controlled through a sampling procedure that balances
        structural compatibility with sequence diversity. The number of generated
        sequences can be adjusted depending on whether the user wants a small set
        of highly confident candidates or a broader exploration of sequence space.

        The imprint ratio parameter regulates how strongly the prediction remains
        biased toward prior sequence information. Lower variability generally
        produces sequences closer to the original scaffold, while higher variability
        encourages exploration of alternative compatible solutions. Biologically,
        this parameter controls the balance between conservation and innovation.

        Two main sampling modes are available. The maximum-confidence strategy
        prioritizes the most probable amino acids at each position, typically
        yielding conservative and stable predictions with lower diversity. The
        probabilistic sampling strategy instead introduces stochastic exploration,
        generating a wider variety of sequences while still maintaining structural
        plausibility. This second mode is often preferred in protein engineering
        campaigns seeking novel sequence solutions.

        Chain and Residue Control

        Advanced options allow selective control over which structural regions are
        redesigned. Entire chains can be excluded from sampling and used instead
        as fixed contextual information during prediction. This capability is
        especially important for multimeric systems where only one partner should
        be redesigned while preserving the remaining assembly.

        The protocol also supports excluding specific amino acids from the sampling
        process. This can be useful for avoiding chemically undesirable residues,
        reducing aggregation propensity, eliminating reactive side chains, or
        adapting sequences to experimental constraints such as synthesis or
        expression limitations.

        Structural Components and imprintRationmental Information

        By default, the protocol incorporates additional structural components
        such as ligands, ions, lipids, and water molecules during prediction.
        These elements may contribute biologically meaningful imprintRationmental
        information that influences sequence compatibility.

        However, users may optionally ignore heteroatoms or solvent molecules.
        This is useful when the surrounding imprintRationment is uncertain, when the
        structure contains crystallographic artifacts, or when the goal is to
        focus exclusively on intrinsic protein features.

        GPU Acceleration and Computational Execution

        The protocol supports both CPU and GPU execution modes. GPU acceleration
        is strongly recommended for large systems or high-throughput sequence
        generation because sampling and downstream structural evaluation can
        become computationally intensive.

        During execution, the protocol automatically prepares the structural
        scaffold, launches the CARBonAra sequence generation engine, collects
        the generated sequences, and organizes the outputs into standardized
        sequence datasets suitable for downstream structural biology workflows.

        Sequence Alignment and Comparative Analysis

        After sequence generation, the protocol automatically aligns all predicted
        sequences together with the original scaffold sequence using Clustal Omega.
        This alignment provides an immediate overview of sequence conservation,
        variability, and consensus regions across the generated ensemble.

        From a biological perspective, the alignment is highly informative because
        conserved positions often correspond to structurally or functionally
        critical residues, while variable positions may indicate regions tolerant
        to mutation or diversification.

        Duplicate sequences are automatically filtered so that downstream analyses
        focus only on unique sequence solutions.

        AlphaFold Structural Evaluation

        Optionally, the protocol can evaluate generated sequences using AlphaFold.
        This step predicts structural confidence metrics for the sampled sequences
        and provides additional validation regarding whether the generated variants
        are likely to maintain the desired fold.

        In multimeric systems, the protocol supports two complementary analyses.
        One evaluates the entire complex, including interface quality and predicted
        interaction confidence. The other focuses specifically on the redesigned
        binder or target chain. This distinction is particularly useful in
        interface engineering and therapeutic binder design.

        Reported metrics may include pLDDT, pTM, and ipTM scores, which together
        provide estimates of local structural confidence, global fold quality,
        and interface reliability. These values help prioritize the most promising
        candidates for experimental validation.

        Output Data and Interpretation

        The protocol produces several categories of outputs. The primary results
        are the generated protein sequences, organized individually and as a
        complete sequence collection. Alignment files summarizing sequence
        conservation are also generated automatically.

        When AlphaFold evaluation is enabled, predicted structural models are
        additionally produced for both complete complexes and redesigned binders.
        These structures can be directly inspected, visualized, or incorporated
        into downstream modeling workflows.

        A summary CSV file is also generated, integrating sequence information,
        CARBonAra confidence scores, and AlphaFold-derived quality metrics. This
        table serves as a practical ranking system for identifying the most
        biologically promising candidates.

        Practical Recommendations

        In most biological applications, it is advisable to begin with moderate
        sequence diversity and a relatively conservative imprint ratio. This
        typically yields structurally reliable candidates while still exploring
        useful sequence variation.

        For binder optimization or interface redesign, preserving contextual
        chains while redesigning only the interacting partner often provides
        more realistic predictions. In contrast, exploratory protein design
        studies may benefit from broader sampling diversity.

        AlphaFold validation is highly recommended when computational resources
        permit, especially before experimental synthesis or screening. Sequences
        with strong structural confidence and favorable interface metrics are
        generally the best candidates for downstream validation.

        Final Perspective

        The Carbonara Sampling Sequence protocol provides a structure-guided
        framework for generating and evaluating alternative protein sequences
        compatible with a target fold. By combining probabilistic sampling,
        structural context, alignment analysis, and optional AlphaFold validation,
        the protocol enables biologically meaningful exploration of protein
        sequence space while maintaining compatibility with experimentally
        relevant structural architectures.
    """
    _label = 'sampling sequence'
    _program = ""
    _version = ""

# -------------------------- DEFINE param functions ----------------------
    METHOD_OPTIONS = ['max', 'sampled']

    def __init__(self, **kwargs):
        """Initialise the protocol and reset complex/binder state flags.

        COMPLEX and BINDER are set during preRequisitesStep based on
        the number of chains in the input structure:
          - COMPLEX = True when the structure has >1 chain
          - BINDER  = True when exactly one chain will be redesigned
                       (single-chain input, or all-but-one chains excluded)
        rows accumulates sequence scores for the final CSV output.
        """
        EMProtocol.__init__(self, **kwargs)
        self.COMPLEX = False
        self.BINDER = False
        self.rows = []

    def _defineParams(self, form):
        """Define the input parameters shown in the protocol form GUI.

        Parameters are grouped into sections:
          - Input: atom structure, number of sequences, imprint ratio
          - Sampling: method (max vs sampled)
          - Chain control: exclude specific chains from sampling
          - Amino acid control: ignore specific amino acids
          - HETATM/water: toggle inclusion of ligands, ions, water
          - GPU: toggle GPU/CPU execution and select GPU ID
        """
        form.addSection(label=Message.LABEL_INPUT)

        # --- Input parameters ---------------------------------------------
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

        form.addParam('imprintRatio', params.FloatParam,
                      default=0.5,
                     validators=[params.GE(0), params.LE(1)],
                      label='Ratio of prior information', important=True,
                       help='Ratio of sequence imprint for sampling'
                            'Value between 0.0 and 1.0. Ratio of the prediction to be used as\n'
                            ' prior information. To get higher variability in the sequence you\n'
                            ' have to choose values close to 1.0. All the positions will contain\n'
                            ' prior information to bias the prediction')

        form.addParam('computeAlphaFold', params.BooleanParam,
                      label="Compute Alphafold scores?",
                      default=False,
                      help="Select 'Yes' if you want to compute AlphaFold scores for sequence\n"
                            "predictions. Take into account that this process cmay take a while.")

        # --- Sampling method ----------------------------------------------
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

        # --- Chain control ------------------------------------------------
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
        form.addParam('selectKnownResidues', params.BooleanParam,
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

        form.addParam('selectUnknownResidues', params.BooleanParam,
                      expertLevel=LEVEL_ADVANCED,
                      label="Select specific RESIDUES for sampling?",
                      default=False,
                      help="Select 'Yes' if you want to select for sampling any specific residues "
                            "(one or more) of the atom structure.\n")

        form.addParam('selectUnknownStructureResidues', params.StringParam,
                      expertLevel=LEVEL_ADVANCED,
                      condition=('selectUnknownResidues==True'),
                      default=None, important=True,
                      label='Included residues',
                      help='"Use the wizard on the right to select the list of resides that will\n'
                           ' be overwritten to the two previous params and re-designed. New\n'
                           ' sequences will be only predicted for these selected positions and\n'
                           ' the rest of the sequence will be used as prior information for the\n'
                           ' prediction."')  
                     # Associate wizard to this parm  
        '''

        # --- Amino acid control -------------------------------------------
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

        # --- HETATM and water control -------------------------------------
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

        # --- GPU execution ------------------------------------------------
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
        """Insert the three processing steps into Scipion's execution queue.

        Step order:
          1. preRequisitesStep (CPU)  - convert CIF, parse chains
          2. processStep (GPU)        - run CARBonAra, alignment, AlphaFold
          3. createOutputStep (CPU)   - register outputs
        """
        self._insertFunctionStep(self.preRequisitesStep, needsGPU=False)
        self._insertFunctionStep(self.processStep, needsGPU=True)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------------------- STEPS functions ------------------------------
    def preRequisitesStep(self):
        """Prepare the input structure before running CARBonAra.

        - Resolves the absolute path of the input atom structure
        - Converts mmCIF to PDB format if needed (CARBonAra requires PDB)
        - Parses chain information to determine if the structure is a
          complex (multiple chains) or a single binder
        - Sets self.COMPLEX and self.BINDER flags for downstream logic
        """
        # get atom structure PDBx/mmCIF file name
        print("self.atomStruct.get(): ", self.atomStruct.get())
        fileName = self.atomStruct.get().getFileName()
        self.atomStructName = os.path.abspath(fileName)
        print("self.atomStructName: ", self.atomStructName)

        # convert CIF to PDB file if needed
        if self.atomStructName.endswith('.cif'):
            atomStructReName = os.path.abspath(
                self._getExtraPath(os.path.basename(
                    os.path.splitext(self.atomStructName)[0] + ".pdb")))
            fromCIFToPDB(
                self.atomStructName, atomStructReName, log=self._log)
            self.atomStructName = atomStructReName

        # get a dict of chains with residues
        h = AtomicStructHandler()
        h.read(self.atomStructName)
        h.getStructure()
        listOfChains, listOfResidues = h.getModelsChains()

        # Extract chain IDs from the first model
        self.chains = list(next(iter(listOfChains.values())))

        # Determine complex vs binder status based on chain count
        if len(self.chains) > 1:
            self.COMPLEX = True
            # If user excluded chains, check if exactly one remains (the binder)
            if self.selectChains and self.selectStructureChains.get():
                selected = [s.strip() for s in self.selectStructureChains.get().split(",")]
                self.binder = [c for c in self.chains if c not in selected]

                if len(self.binder) == 1:
                    self.BINDER = True
                    # Store the index of the single redesigned chain
                    self.idx = self.chains.index(str(self.binder[0]))
        elif len(self.chains) == 1:
            # Single-chain input is always treated as a binder
            self.BINDER = True
            self.idx = 0
        self._log.info(f"self.COMPLEX: {self.COMPLEX}")
        self._log.info(f"self.BINDER: {self.BINDER}")

    def processStep(self):
        """Main processing step: build CLI args, run CARBonAra, align sequences.

        This step:
          1. Builds the command-line arguments for the carbonara binary
          2. Runs CARBonAra via self.runJob()
          3. Performs ClustalO alignment of generated sequences
          4. Optionally runs ColabFold/AlphaFold for structural scoring
          5. Extracts and merges scores into a summary CSV
        """
        args = []

        # --- CARBonAra CLI arguments --------------------------------------

        # Number of sequences to generate
        args.extend(["--num_sequences", str(self.numSamples.get())])

        # ratio of prior information (0.0-1.0, higher = more variability)
        args.extend(["--imprint_ratio", str(self.imprintRatio.get())])

        # sampling method: 'max' (high confidence) or 'sampled' (higher diversity)
        args.extend(["--sampling_method", self.METHOD_OPTIONS[self.bSampled.get()]])

        # excluded chains (chains used as fixed context, not redesigned)
        if self.selectChains and \
           self.selectStructureChains.get():
            chains = str(self.selectStructureChains.get())
            args.extend(["--known_chains", chains.replace(" ", "")])
        '''
        TODO: Restore this parameter when it works in carbonara method

        # excluded residues
        if (self.selectKnoumResidues==True and
            self.selectKnownStructureResidues.get()):
            knownResidues = str(self.selectKnownResidues.get())
            args.extend(["--known_positions", knownResidues])

        # included residues  
        if (self.selectUnknounResidues==True and 
            self.selectUnknownStructureResidues.get()):
            unknownResidues = str(self.selectUnknownStructureResidues.get())
            args.extend(["--unknown_positions", unknownResidues]) 
        '''

        # ignored amino acids (excluded from sampling entirely)
        if self.ignoreAminoacids and\
           self.selectIgnoredAminoacids.get():
            ignoredResidues = str(self.selectIgnoredAminoacids.get())
            args.extend(["--ignored_amino_acids",
                         ignoredResidues.replace(" ", "")])

        # HETATM toggle (ligands, ions, lipids)
        if self.ignoreHetatm:
            args.extend(["--ignore_hetatm"])

        # water toggle
        if self.ignoreWater:
            args.extend(["--ignore_water"])

        # device selection: CPU or GPU (with CUDA_VISIBLE_DEVICES)
        if not self.useGpu:
            args.extend(["--device", "cpu"])
            self._log.info("I'm using CPU ")
        else:
            args.extend(["--device", "cuda"])
            os.environ["CUDA_VISIBLE_DEVICES"] = ("%s" % self.getGpuList()[0])
            self._log.info("I'm using GPU " + ("%s" % self.getGpuList()[0]))

        # input PDB file path (positional argument)
        args.extend([self.atomStructName])

        # output directory: <extra_path>/carbonara_results
        self.dir_path = os.path.abspath(self._getExtraPath())
        subdirectory = "carbonara_results"
        self.subdir_path = os.path.join(self.dir_path, subdirectory)
        args.extend([self.subdir_path])

        # Join args into a single string for shell execution
        args_str = ' '.join(args)

        # Run CARBonAra via the plugin's conda-activated command
        self.runJob(Plugin.getCarbonaraCmd(), args_str)

        # --- Post-processing: alignment and scoring -----------------------

        # Align all generated sequences with the original using ClustalO
        self.clustalOAlignment()

        # Generate a filtered alignment summary showing only consensus symbols
        summ_align_file = os.path.abspath(
            self._getExtraPath("clustal_summary.aln"))
        self.filter_file_lines(self.outFile, summ_align_file)

        # --- Optional AlphaFold/ColabFold scoring -------------------------
        if self.computeAlphaFold:
            # AlphaFold prediction for complex (two or more chains)
            # Produces pLDDT, pTM, ipTM scores (ipTM is null for single-chain)
            if self.COMPLEX:
                self._log.info("Running colabfold_batch, this may take a while...")
                merged_multimer_folder_path = os.path.abspath(
                    self._getExtraPath("merged.fasta"))
                alphafold_subdirectory_multimer =\
                    "alphafold_predictions_multimer"
                self.multimer_folder_path = os.path.join(
                    self.dir_path, alphafold_subdirectory_multimer)
                self.compute_alphafold(
                    merged_multimer_folder_path, self.multimer_folder_path)

            # AlphaFold prediction for the redesigned binder chain only
            # (or for single-chain inputs)
            # Produces pLDDT, pTM scores
            if self.BINDER:
                alphafold_subdirectory_binder = "alphafold_predictions_binder"
                self.binder_folder_path = os.path.join(
                    self.dir_path, alphafold_subdirectory_binder)
                if len(self.chains) > 1:
                    # Complex case: extract the binder chain from the merged FASTA
                    merged_binder_file_path = os.path.abspath(
                        self._getExtraPath("merged_binder.fasta"))
                    self.extract_binder_chain_from_aln_file(
                        merged_multimer_folder_path,
                        self.idx,
                        merged_binder_file_path)
                    self.compute_alphafold(
                        merged_binder_file_path, self.binder_folder_path)
                if len(self.chains) == 1:
                    # Single-chain case: use the full merged FASTA directly
                    merged_multimer_folder_path = os.path.abspath(
                        self._getExtraPath("merged.fasta"))
                    self.compute_alphafold(
                        merged_multimer_folder_path, self.binder_folder_path)

        # --- Score extraction and CSV generation --------------------------

        # Extract CARBonAra scores from output FASTA files
        output_scores_csv_file = os.path.join(self.dir_path, "sorted_scores.csv")
        self.extract_scores_from_fasta(
            self.subdir_path, output_scores_csv_file)

        # If AlphaFold was run, merge its metrics into the CSV
        if self.computeAlphaFold:
            metrics_complex = None
            metrics_binder = None

            # Extract complex metrics from colabfold log if available
            if (self.COMPLEX and
               os.path.exists(os.path.join(self.multimer_folder_path, "log.txt"))):
               log_file = os.path.join(self.multimer_folder_path, "log.txt")
               metrics_complex= self.extract_metrics_from_log_complex(log_file)

            # Extract binder metrics from colabfold log if available
            if (self.BINDER and
                os.path.exists(os.path.join(self.binder_folder_path, "log.txt"))):
                log_file = os.path.join(self.binder_folder_path, "log.txt")
                metrics_binder = self.extract_metrics_from_log_binder(log_file)

            # Merge CARBonAra scores with AlphaFold metrics and write final CSV
            self.merge_and_write_csv(
                self.rows, output_scores_csv_file,
                    metrics_complex=metrics_complex,
                    metrics_binder=metrics_binder)

    def createOutputStep(self):
        """Register generated sequences and optional AlphaFold structures as Scipion outputs.

        Outputs:
          - Individual Sequence objects for each unique generated sequence
          - A SetOfSequences containing all unique sequences
          - Optional AtomStruct objects for AlphaFold-predicted structures
            (complex and/or binder, depending on input)
        """
        outputs = {}

        # Create the SetOfSequences container
        setSeq = SetOfSequences()
        outputSequences = setSeq.create(outputPath=self._getPath())

        sequences = []
        # Get list of unique sequence file paths (deduplicated by sequence content)
        filepath_list = self.sel_files_seq_unique(self.subdir_path)

        # --- Register individual sequences --------------------------------
        # Iterate over all .fasta files in carbonara_results/
        for filename in sorted(os.listdir(self.subdir_path)):
            if filename.endswith(".fasta"):
                path = os.path.join(self.subdir_path, filename)
                sequence = self.extract_sequence_from_file(path)
                label = os.path.splitext(filename)[0]

                # Only register unique sequences (skip duplicates)
                if path in filepath_list:
                    seq = Sequence()
                    seq.setSeqName(path)
                    seq.setId(label)
                    seq.setSequence(sequence)
                    seq.setObjLabel(label)
                    # Register as individual output with filename as key
                    keyword = filename.split(".fasta")[0]
                    kwargs = {keyword: seq}
                    self._defineOutputs(**kwargs)

                    # Add to the set
                    sequences.append(seq)

        # Append all sequences to the SetOfSequences
        for seq in sequences:
            outputSequences.append(seq)

        self._defineOutputs(outputSequences=outputSequences)

        # total number of unique sequences
        num_sequences = len(sequences)
        self._log.info(f"Total number of unique sequences: {num_sequences}")

        # --- Register AlphaFold-predicted structures ----------------------
        if self.computeAlphaFold:
            # Register complex structures (multi-chain predictions)
            if self.COMPLEX:
                for filename in sorted(os.listdir(self.multimer_folder_path)):
                    if filename.endswith(".pdb") or filename.endswith(".cif"):
                        path = self.multimer_folder_path / filename
                        pdb = AtomStruct()
                        pdb.setFileName(path)
                        keyword = filename.split("_unrelaxed_")[0] + "_complex"
                        outputs[keyword] = pdb

            # Register binder structures (single-chain predictions)
            if self.BINDER:
                for filename in sorted(os.listdir(self.binder_folder_path)):
                    if filename.endswith(".pdb") or filename.endswith(".cif"):
                        path = os.path.join(self.binder_folder_path, filename)
                        pdb = AtomStruct()
                        pdb.setFileName(path)
                        keyword = filename.split("_unrelaxed_")[0] + "_binder"
                        outputs[keyword] = pdb

            if outputs:
                self._defineOutputs(**outputs)



    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        """Validate protocol configuration before execution.

        Checks:
          - Only one GPU is selected (CARBonAra uses a single GPU)
          - clustalo binary is available in the carbonara conda env
        Returns a list of error messages (empty if validation passes).
        """
        errors = []

        gpus = self.getGpuList()

        if len(gpus) > 1:
            errors.append('Only one GPU can be used.')

        # Check that CLUSTALO program exists in the carbonara conda env
        result = subprocess.run(
            ["conda", "env", "list", "--json"],
             capture_output=True,
             text=True,
             check=True
        )
        envs = json.loads(result.stdout)["envs"]
        env_path = next(p for p in envs if p.endswith(f"/{conda_env}"))
        bin_dir = env_path + ("/Scripts" if os.name == "nt" else "/bin")
        self._log.info(f"Checking for Clustal Omega in: {bin_dir}")
        if not (self.is_tool(os.path.join(bin_dir,CLUSTALO))):
            errors.append("Clustal-omega program missing.\n "
                          "You need it to run this program.\n"
                          "Please install Clustal-omega:\n"
                          "     sudo apt-get install clustalo\n")

        return errors

    def _summary(self):
        """Generate a text summary of protocol results for display in the GUI.

        Includes:
          - Total number of predicted sequences
          - Number of unique sequences (after deduplication)
          - Path to the ClustalO alignment file
          - The ClustalO alignment summary (consensus symbols only)
        """
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
        """Return the list of citation keys for this protocol."""
        return ['Krapp2024-dw']

    def clustalOAlignment(self):
        """Perform multiple sequence alignment using Clustal Omega.

        Steps:
          1. Merge the original sequence and all sampled sequences into a
             single FASTA file (merged.fasta)
          2. Run ClustalO to generate the alignment (clustal.aln)
        The alignment is later filtered to produce a summary with consensus
        symbols only.
        """
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
        """Merge the original sequence and all sampled sequences into one FASTA file.

        The merged file (merged.fasta) is used as input for ClustalO alignment.
        Duplicate sequences are excluded via sel_files_seq_unique().

        Returns the path to the merged FASTA file.
        """
        inFile =  os.path.abspath(self._getExtraPath("merged.fasta"))
        # Extract the original protein sequence(s) from the input structure
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
        """Extract protein sequences from a PDB/mmCIF file using Biopython.

        Parses the first model of the structure, extracts all standard amino
        acid residues from each chain, and returns them as a colon-separated
        string (e.g. 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH').

        Multiple chains are joined with ':' to preserve chain boundaries.
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", pdb_file)

        chain_sequences = []

        for model in structure:
            for chain in model:
                # Filter to standard amino acids only
                residues = [res for res in chain if is_aa(res, standard=True)]
                if residues:
                    sequence = ''.join(seq1(res.get_resname()) for res in residues)
                    chain_sequences.append(sequence)
            break  # Only process the first model

        return ':'.join(chain_sequences)

    def extract_scores_from_fasta(self, folder_path, output_scores_csv_file):
        """Extract CARBonAra scores from output FASTA files and write a summary CSV.

        Parses the FASTA headers to extract the 'score=...' value, deduplicates
        sequences, and writes a CSV with columns: file_name, sequence, score.

        If AlphaFold is not enabled, the CSV is sorted by score (highest first).
        If AlphaFold is enabled, the CSV is written later by merge_and_write_csv()
        after merging with AlphaFold metrics.

        Populates self.rows with [file_name, sequence, score] entries.
        """
        # include starting sequence as first row with score=1 (reference)
        name = os.path.basename(
            os.path.splitext(self.atomStructName)[0])
        self.rows.append([name, self.wholeSequence, 1])

        #include all the sampled sequences without duplicates
        filepath_list = self.sel_files_seq_unique(folder_path)

        for filename in sorted(os.listdir(folder_path)):
            if filename.endswith(".fasta"):
                filepath = os.path.join(folder_path, filename)

                if filepath in filepath_list:
                    seq_name = os.path.splitext(filename)[0]
                    for record in SeqIO.parse(filepath, "fasta"):
                        seq_str = str(record.seq)

                        #extract score from FASTA header (format: "score=X.XX")
                        header = record.description
                        score = None
                        for part in header.split(","):
                            if "score" in part:
                                score = float(part.split("=")[1].strip())
                                self.rows.append([seq_name, seq_str, score])
                                break

        if not self.computeAlphaFold:
            # sort rows by score descending and write CSV immediately
            self.rows.sort(key=lambda x: x[2], reverse=True)

           # write to CSV
            with open(output_scores_csv_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(["file_name", "sequence", "CARBonAra score"]) # header row
                # then all extracted rows
                writer.writerows(self.rows)

            self._log.info(f"CSV file generated: {output_scores_csv_file} with {len(self.rows)} unique sequences")

    def filter_file_lines(self, filepath, newfilepath):
        """Filter a ClustalO alignment file to keep only the reference sequence
        and consensus symbols.

        The filtered output (clustal_summary.aln) shows:
          - The original/reference sequence (first_label)
          - Consensus lines (* for identical, : for conserved, . for semi-conserved)
          - CLUSTAL header and empty separator lines

        This produces a compact view of alignment quality without cluttering
        the display with all sampled sequences.
        """
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

                    # Keep consensus lines (start with *, :, or .)
                    if stripped[0] in "*:." or line.startswith(" "):
                        filtered_lines.append(line)
                        continue

                    # Process sequence lines (format: "label  sequence")
                    parts = line.rstrip().split()
                    if len(parts) == 2:
                        label, seq = parts
                        if first_label is None:
                            # First sequence encountered is the reference
                            first_label = label
                            filtered_lines.append(line)
                        elif first_label and label == first_label:
                            # Keep all blocks of the reference sequence
                            filtered_lines.append(line)
                    # Skip other sequence labels

            # New file with valid lines
            with open(newfilepath, "w") as f:
                f.writelines(filtered_lines)

        except Exception as e:
            self._log.error(f"Error processing file: {e}")

    def extract_sequence_from_file(self, filepath):
        """Extract the raw sequence string from a FASTA file.

        Skips header lines (starting with '>') and concatenates all
        sequence lines into a single string.
        """
        sequence_lines = []
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('>'):
                    continue
                sequence_lines.append(line)
        return ''.join(sequence_lines)

    def sel_files_seq_unique(self, folder_path):
        """Return a list of FASTA file paths with unique sequences.

        Iterates over all .fasta files in folder_path, extracts the sequence
        from each, and returns only the paths of files whose sequences have
        not been seen before. This deduplicates sequences that CARBonAra
        may generate multiple times.
        """
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

    def compute_alphafold(self, in_folder_path, out_folder_path):
        """Run ColabFold/AlphaFold predictions on sequences in a FASTA file.

        Executes colabfold_batch in the colabfold conda environment to predict
        structures for all sequences in in_folder_path, writing results to
        out_folder_path.

        Raises an exception if colabfold_batch exits with a non-zero code.
        """
        command = [
            "conda", "run", "-n", "colabfold",
            "colabfold_batch", "--num-models", "1",
            in_folder_path, out_folder_path
        ]

        self._log.info(greenStr(f"Running command: {' '.join(command)}"))
        # self.runJob(command)
        result = subprocess.run(
            command,
            capture_output=True,
            text=True
        )

        stdout = result.stdout
        stderr = result.stderr
        return_code = result.returncode
        self._log.info(f"Command stdout: {stdout}")
        self._log.info(f"Command stderr: {stderr}")
        self._log.info(f"Command return code: {return_code}")

        if return_code != 0:
            raise Exception(
                f"colabfold_batch failed (exit code {return_code}) for "
                f"input '{in_folder_path}'. stderr: {stderr}")

    def extract_binder_chain_from_aln_file(self, in_folder_path, chain_index, out_folder_path):
        """Extract a specific chain from multi-chain FASTA sequences.

        Input FASTA files contain sequences with chains separated by ':'
        (e.g. 'MVLSPADKTN:GSAQVKGHGK'). This method extracts the chain at
        chain_index (0-based) and writes it to a new FASTA file with an
        updated header indicating the chain ID.

        Parameters:
          - in_folder_path: path to the input FASTA file (multi-chain)
          - chain_index: 0-based index of the chain to extract (0=A, 1=B, ...)
          - out_folder_path: path to the output FASTA file (single chain)
        """
        with open(in_folder_path) as fin, open(out_folder_path, "w") as fout:
            header = None
            for line in fin:
                line = line.strip()
                if line.startswith(">"):
                    # Save the header line (sequence identifier)
                    header = line
                else:
                    # Split the sequence into chains using ':'
                    chains = line.split(":")
                    if chain_index < len(chains):
                        seq = chains[chain_index]
                        # Write the new header with chain index
                        fout.write(f"{header}_chain{self.chains[chain_index]}\n")
                        # Write the extracted sequence
                        fout.write(seq + "\n")

    def extract_metrics_from_log_complex(self, log_file):
        """Extract AlphaFold metrics for complex (multi-chain) predictions.

        Parses the colabfold_batch log file to extract pLDDT, pTM, and ipTM
        scores for each sequence. ipTM (interface pTM) is only meaningful for
        multi-chain complexes.

        Returns a list of [seq_name, pLDDT, pTM, ipTM] entries.
        """
        results = []

        with open(log_file, "r") as f:
            lines = f.readlines()

        seq_name = None
        for line in lines:

            # Detect sequence name from "Query N/M: seq_name" lines
            if "Query" in line:
                match = re.search(r"Query\s+\d+/\d+:\s+(\S+)", line)
                if match:
                    seq_name = match.group(1)

            # Extract scores from "rank_001" lines containing pLDDT/pTM/ipTM
            if "rank_001" in line and seq_name:
                plddt = re.search(r"pLDDT=([\d\.]+)", line)
                ptm   = re.search(r"pTM=([\d\.]+)", line)
                iptm  = re.search(r"ipTM=([\d\.]+)", line)

                results.append([
                    seq_name,
                    float(plddt.group(1)) if plddt else None,
                    float(ptm.group(1)) if ptm else None,
                    float(iptm.group(1)) if iptm else None
                ])

                # Reinitiate seq_name to avoid association to other line by error
                seq_name = None

        return results

    def extract_metrics_from_log_binder(self, log_file):
        """Extract AlphaFold metrics for binder (single-chain) predictions.

        Parses the colabfold_batch log file to extract pLDDT and pTM scores
        for each sequence. Does not extract ipTM since it's only meaningful
        for multi-chain complexes.

        Returns a list of [seq_name, pLDDT, pTM] entries.
        """
        results = []

        with open(log_file, "r") as f:
            lines = f.readlines()

        seq_name = None
        for line in lines:

            # Detect sequence name from "Query N/M: seq_name" lines
            if "Query" in line:
                match = re.search(r"Query\s+\d+/\d+:\s+(\S+)", line)
                if match:
                    seq_name = match.group(1)

            # Extract scores from "rank_001" lines containing pLDDT/pTM
            if "rank_001" in line and seq_name:
                plddt = re.search(r"pLDDT=([\d\.]+)", line)
                ptm   = re.search(r"pTM=([\d\.]+)", line)

                results.append([
                    seq_name,
                    float(plddt.group(1)) if plddt else None,
                    float(ptm.group(1)) if ptm else None,
                ])

                # Reinitiate seq_name to avoid association to other line by error
                seq_name = None

        return results


    def merge_and_write_csv(self, seq_scores, output_file,
                            metrics_complex=None, metrics_binder=None):
        """Merge CARBonAra scores with AlphaFold metrics and write the final CSV.

        Combines three data sources into a single ranked CSV:
          - seq_scores: CARBonAra scores (file_name, sequence, score)
          - metrics_complex: AlphaFold metrics for complex predictions
            (file_name, pLDDT, pTM, ipTM)
          - metrics_binder: AlphaFold metrics for binder predictions
            (file_name_chain, pLDDT, pTM)

        Sorting priority (highest to lowest):
          - If both complex and binder metrics exist: ipTM_complex, pTM_complex, pTM_binder
          - If only complex metrics: ipTM_complex, pTM_complex
          - If only binder metrics: pTM_binder
          - Otherwise: CARBonAra score

        Writes CSV with columns: file_name, sequence, score, pLDDT_complex,
        pTM_complex, ipTM_complex, pLDDT_binder, pTM_binder
        """
        # Convert seq_scores to dict keyed by file_name for easy lookup
        merged = {}

        for row in seq_scores:
            # self._log.info(f"row: {row}")
            file_name, sequence, score = row
            merged[file_name] = {"file_name": file_name,
                                 "sequence": sequence,
                                 "score": score}

        # Merge complex metrics (pLDDT, pTM, ipTM) into the dict
        if metrics_complex:
            for row in metrics_complex:
                file_name, plddt, ptm, iptm = row

                if file_name in merged:
                    merged[file_name]["pLDDT_complex"] = plddt
                    merged[file_name]["pTM_complex"] = ptm
                    merged[file_name]["ipTM_complex"] = iptm

        # Merge binder metrics (pLDDT, pTM) into the dict
        # Binder filenames have "_chainX" suffix that must be stripped
        if metrics_binder:
            for row in metrics_binder:

                file_name_chain, plddt, ptm = row
                try:
                    file_name = file_name_chain.split("_chain")[0]
                except:
                    file_name = file_name_chain

                if file_name in merged:
                    merged[file_name]["pLDDT_binder"] = plddt
                    merged[file_name]["pTM_binder"] = ptm

        # convert to list of rows for sorting
        rows = list(merged.values())

        # Determine sorting strategy based on available metrics
        has_complex = any("ipTM_complex" in r for r in rows)
        has_binder  = any("pTM_binder" in r for r in rows)

        if has_complex and has_binder:
            # Sort by ipTM_complex (interface quality), then pTM_complex, then pTM_binder
            rows.sort(key=lambda r: (
                r.get("ipTM_complex", -1),
                r.get("pTM_complex", -1),
                r.get("pTM_binder", -1)
            ), reverse=True)
            sort_msg = "orderd by ipTM_complex, then pTM_complex and finally pTM_binder"
        elif has_complex:
            # Sort by ipTM_complex, then pTM_complex
            rows.sort(key=lambda r: (
                r.get("ipTM_complex", -1),
                r.get("pTM_complex", -1)
            ), reverse=True)
            sort_msg = "ordered by ipTM_complex and then pTM_complex"
        elif has_binder:
            # Sort by pTM_binder only
            rows.sort(key=lambda r: r.get("pTM_binder", -1), reverse=True)
            sort_msg = "ordered by pTM_binder"
        else:
            # Fall back to CARBonAra score
            rows.sort(key=lambda r: r["score"], reverse=True)
            sort_msg = "ordered by score"

        # determine header (all possible columns; missing values will be None)
        header = [
            "file_name", "sequence", "score",
            "pLDDT_complex", "pTM_complex", "ipTM_complex",
            "pLDDT_binder", "pTM_binder"
        ]

        # write CSV, filling missing columns with None
        with open(output_file, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=header)
            writer.writeheader()
            for r in rows:
                for h in header:
                    if h not in r:
                        r[h] = None
                writer.writerow(r)
        self._log.info(f"CSV generated: {output_file} with {len(merged)} sequences, {sort_msg}.")
