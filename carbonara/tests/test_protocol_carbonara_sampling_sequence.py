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
Integration tests for the CarbonaraSamplingSequence protocol.

Tests run through Scipion's test framework (not pytest) and require:
  - A working Scipion installation with the carbonara plugin
  - The carbonara and colabfold conda environments built
  - The scipion-em-chimera plugin (for ChimeraProtOperate)
  - Network access (imports PDB 5ni1 from RCSB)

Run with: scipion3 tests carbonara.tests.test_protocol_carbonara_sampling_sequence
"""

import os
from ..protocols import CarbonaraSamplingSequence
from chimera.protocols import ChimeraProtOperate
from pyworkflow.tests import BaseTest, setupTestProject
from pwem.protocols.protocol_import import (ProtImportPdb)


class TestImportBase(BaseTest):
    """Base test class that sets up a Scipion test project."""
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)


class TestImportData(TestImportBase):
    """Helper class that provides PDB import functionality for other tests."""
    def _importStructurePDB(self, pdbID="5ni1"):  # Haemoglobin
        """Import a PDB structure from the RCSB by ID.

        Default is 5ni1 (Haemoglobin), a tetrameric protein useful for
        testing multi-chain scenarios.
        """
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': pdbID
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import structure Haemoglobin\n %s' % pdbID)
        self.launchProtocol(protImportPDB)
        self.assertTrue(protImportPDB.outputPdb.getFileName())
        return protImportPDB.outputPdb


class TestCarbonaraSamplingSequence(TestImportData):
    """Integration tests for the CarbonaraSamplingSequence protocol.

    Tests cover:
      - Default parameters with multi-chain input (testCarbonara1)
      - Chain exclusion scenarios
      - AlphaFold scoring integration
      - Single-chain input via Chimera chain extraction (testCarbonara2)
      - Various imprint ratios and sampling methods
      - Amino acid exclusion
      - CPU vs GPU execution
    """

    def testCarbonara1(self):
        """Test CARBonAra with multi-chain input (Haemoglobin, 4 chains).

        Scenarios:
          1. Default parameters, 3 sequences
          2. Exclude chain B from sampling
          3. Exclude chain B + AlphaFold scoring
          4. Exclude chains B, D (known CARBonAra bug)
          5. Exclude chains B, D + AlphaFold (known bug)
          6. Exclude chains B, C, D + AlphaFold (known bug)
        """
        # Scenario 1: Default parameters, 3 sequences
        # This test checks that carbonara runs with default parameters
        #     except number of sequences (3)
        #
        print("Run Carbonara Sampling Sequence protocol from imported pdb file " \
              "with default parameters and 3 predicted sequences")

        # import PDB
        Haemoglobin_PDB = self._importStructurePDB()
        self.assertTrue(os.path.exists(Haemoglobin_PDB.getFileName()))

        args = {'atomStruct': Haemoglobin_PDB,
                'numSamples': 3
                }

        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara default params')
        self.launchProtocol(protCarbonaraSamplingSequence)
        outPuts = ['5ni1_0', '5ni1_1', '5ni1_2']
        t = protCarbonaraSamplingSequence.iterOutputAttributes()
        for output, ProtOutput in zip(outPuts, 
                                      protCarbonaraSamplingSequence.
                                      iterOutputAttributes()):
            self.assertEqual(str(ProtOutput[0]), output)
            print(str(ProtOutput[0]))


        # Scenario 2: Exclude chain B from sampling
        # This test checks that carbonara runs excluding a chain from sampling
        #
        print("Run Carbonara Sampling Sequence protocol from imported pdb file " \
              "and excluding chain B from sampling")

        args = {'atomStruct': Haemoglobin_PDB,
                'numSamples': 3,
                'selectChains' : True,
                'selectStructureChains' : 'B'
                }

        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara excluded chain B')
        self.launchProtocol(protCarbonaraSamplingSequence)
        outPuts = ['5ni1_0', '5ni1_1', '5ni1_2']
        for output, ProtOutput in zip(outPuts, 
                                      protCarbonaraSamplingSequence.
                                      iterOutputAttributes()):
            self.assertEqual(str(ProtOutput[0]), output)
            print(str(ProtOutput[0]))

        # Scenario 3: Exclude chain B + AlphaFold scoring
        # This test checks that carbonara runs excluding a chain from sampling
        # and launches AlphaFold
        #
        print("Run Carbonara Sampling Sequence protocol from imported pdb file " \
              "and excluding chain B from sampling and launches AlphaFold")

        args = {'atomStruct': Haemoglobin_PDB,
                'numSamples': 3,
                'computeAlphaFold' : True,
                'selectChains' : True,
                'selectStructureChains' : 'B'
                }

        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara excluded chain B alpahfold')
        self.launchProtocol(protCarbonaraSamplingSequence)
        outPuts = ['5ni1_0', '5ni1_1', '5ni1_2']
        for output, ProtOutput in zip(outPuts, 
                                      protCarbonaraSamplingSequence.
                                      iterOutputAttributes()):
            self.assertEqual(str(ProtOutput[0]), output)
            print(str(ProtOutput[0]))


        # Scenario 4: Exclude chains B, D (known CARBonAra bug)
        # This test checks that carbonara runs excluding two chains from sampling
        #
        print("Run Carbonara Sampling Sequence protocol from imported pdb file " \
              "and excluding chains B and D from sampling" \
              "We except that this test fails due to a carbonara bugg")


        args = {'atomStruct': Haemoglobin_PDB,
                'numSamples': 3,
                'selectChains' : True,
                'selectStructureChains' : 'B, D',
                }

        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara excluded chains B, D')
        self.launchProtocol(protCarbonaraSamplingSequence)
        outPuts = ['5ni1_0', '5ni1_1', '5ni1_2']
        for output, ProtOutput in zip(outPuts, 
                                      protCarbonaraSamplingSequence.
                                      iterOutputAttributes()):
            self.assertEqual(str(ProtOutput[0]), output)
            print(str(ProtOutput[0]))

        # Scenario 5: Exclude chains B, D + AlphaFold (known bug)
        # This test checks that carbonara runs excluding two chains from sampling
        # and launches AlphaFold
        #
        print("Run Carbonara Sampling Sequence protocol from imported pdb file " \
              "and excluding chains B and D from sampling" \
              " and launches AlphaFold"
              "We except that this test fails due to a carbonara bugg")


        args = {'atomStruct': Haemoglobin_PDB,
                'numSamples': 3,
                'computeAlphaFold' : True,
                'selectChains' : True,
                'selectStructureChains' : 'B, D',
                }

        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara excluded chains B, D, alphafold')
        self.launchProtocol(protCarbonaraSamplingSequence)
        outPuts = ['5ni1_0', '5ni1_1', '5ni1_2']
        for output, ProtOutput in zip(outPuts, 
                                      protCarbonaraSamplingSequence.
                                      iterOutputAttributes()):
            self.assertEqual(str(ProtOutput[0]), output)
            print(str(ProtOutput[0]))

        # Scenario 6: Exclude chains B, C, D + AlphaFold (known bug)
        # This test checks that carbonara runs excluding three chains from sampling
        # and launches AlphaFold
        #
        print("Run Carbonara Sampling Sequence protocol from imported pdb file " \
              "and excluding chains B, C and D from sampling" \
              " and launches AlphaFold"
              "We except that this test fails due to a carbonara bugg")


        args = {'atomStruct': Haemoglobin_PDB,
                'numSamples': 3,
                'computeAlphaFold' : True,
                'selectChains' : True,
                'selectStructureChains' : 'B, C, D',
                }

        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara excluded chains B,C,D, alphafold')
        self.launchProtocol(protCarbonaraSamplingSequence)
        outPuts = ['5ni1_0', '5ni1_1', '5ni1_2']
        for output, ProtOutput in zip(outPuts, 
                                      protCarbonaraSamplingSequence.
                                      iterOutputAttributes()):
            self.assertEqual(str(ProtOutput[0]), output)
            print(str(ProtOutput[0]))

    def testCarbonara2(self):
        """Test CARBonAra with single-chain input extracted via Chimera.

        Uses ChimeraProtOperate to extract chain A from Haemoglobin, then
        runs CARBonAra with various parameter combinations:
          - imprint ratio 0 and 1, sampling method 'sampled'
          - imprint ratio 0 and 1, sampling method 'max'
          - ignored amino acids (Cys, Trp)
          - CPU vs GPU execution
          - AlphaFold scoring integration
        """
        # This test checks that carbonara runs with default parameters
        #     except number of sequences (3) starting from a unique chain of the atom structure
        #
        print("Run Carbonara Sampling Sequence protocol from a chimerax operate protocol " \
              "that select only one chain of the atom structure (chain A) with default "
              "parameters and generates 3 predicted sequences")

        # import PDB
        Haemoglobin_PDB = self._importStructurePDB()
        self.assertTrue(os.path.exists(Haemoglobin_PDB.getFileName()))

        # extract chain A with chimerax-operate protocol
        extraCommands = ""
        extraCommands += "select ~/A \n" 
        extraCommands += "delete sel \n"
        extraCommands += "scipionwrite #2 " \
                         "prefix DONOTSAVESESSION_ \n"
        extraCommands += "exit \n"

        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': Haemoglobin_PDB
                }
        protChimera = self.newProtocol(ChimeraProtOperate, **args)
        protChimera.setObjLabel('chimera operate chain A extracted')
        self.launchProtocol(protChimera)                
        PDB_output = eval("protChimera.DONOTSAVESESSION_Atom_struct__2_%06d"
                          % protChimera.getObjId())
        self.assertTrue(os.path.exists(PDB_output.getFileName()))

        # Sub-test 1: imprint ratio 0, sampling method 'sampled'
        # run carbonara 3 seq, inprint ratio 0, sampling method sampled 
        args = {'atomStruct': PDB_output,
                'numSamples': 3,
                'imprintRatio' : 0
                }
        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara inp_ratio 0 sampled')
        self.launchProtocol(protCarbonaraSamplingSequence)

        outPuts = []
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_0")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_1")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_2")

        ProtOutputs = []
        for t in protCarbonaraSamplingSequence.iterOutputAttributes():
            name = t[0]
            if name.startswith("DONOTSAVESESSION_Atom_struct__2_000042"):
                ProtOutputs.append(name)
        
        for output, ProtOutput in zip(outPuts, ProtOutputs):
            self.assertEqual(output, ProtOutput)

        # Sub-test 2: imprint ratio 1, sampling method 'sampled'
        # run carbonara 3 seq, inprint ratio 1, sampling method sampled 
        args = {'atomStruct': PDB_output,
                'numSamples': 3,
                'imprintRatio' : 1
                }
        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara inp_ratio 1 sampled')
        self.launchProtocol(protCarbonaraSamplingSequence)
        
        outPuts = []
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_0")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_1")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_2")

        ProtOutputs = []
        for t in protCarbonaraSamplingSequence.iterOutputAttributes():
            name = t[0]
            if name.startswith("DONOTSAVESESSION_Atom_struct__2_000042"):
                ProtOutputs.append(name)
        
        for output, ProtOutput in zip(outPuts, ProtOutputs):
            self.assertEqual(output, ProtOutput)

        # Sub-test 3: imprint ratio 0, sampling method 'max'
        # run carbonara 3 seq, inprint ratio 0, sampling method max
        args = {'atomStruct': PDB_output,
                'numSamples': 3,
                'imprintRatio' : 0,
                'bSampled' : 0
                }
        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara inp_ratio 0 max')
        self.launchProtocol(protCarbonaraSamplingSequence)

        outPuts = []
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_0")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_1")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_2")

        ProtOutputs = []
        for t in protCarbonaraSamplingSequence.iterOutputAttributes():
            name = t[0]
            if name.startswith("DONOTSAVESESSION_Atom_struct__2_000042"):
                ProtOutputs.append(name)
        
        for output, ProtOutput in zip(outPuts, ProtOutputs):
            self.assertEqual(output, ProtOutput)

        # Sub-test 4: imprint ratio 1, sampling method 'max'
        # run carbonara 3 seq, inprint ratio 1, sampling method max
        args = {'atomStruct': PDB_output,
                'numSamples': 3,
                'imprintRatio' : 1,
                'bSampled' : 0
                }
        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara inp_ratio 1 max')
        self.launchProtocol(protCarbonaraSamplingSequence)
        
        outPuts = []
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_0")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_1")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_2")

        ProtOutputs = []
        for t in protCarbonaraSamplingSequence.iterOutputAttributes():
            name = t[0]
            if name.startswith("DONOTSAVESESSION_Atom_struct__2_000042"):
                ProtOutputs.append(name)
        
        for output, ProtOutput in zip(outPuts, ProtOutputs):
            self.assertEqual(output, ProtOutput)

        # Sub-test 5: imprint ratio 1, sampling method 'sampled', ignore Cys
        # run carbonara 3 seq, inprint ratio 1, sampling method samlpled, ignored aminoacid Cys
        args = {'atomStruct': PDB_output,
                'numSamples': 3,
                'imprintRatio' : 1,
                'bSampled' : 1,
                'ignoreAminoacids' : True,
                'selectIgnoredAminoacids' : "C"
                }
        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara inp_ratio 1 sampled C ignored')
        self.launchProtocol(protCarbonaraSamplingSequence)
        
        outPuts = []
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_0")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_1")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_2")

        ProtOutputs = []
        for t in protCarbonaraSamplingSequence.iterOutputAttributes():
            name = t[0]
            if name.startswith("DONOTSAVESESSION_Atom_struct__2_000042"):
                ProtOutputs.append(name)
        
        for output, ProtOutput in zip(outPuts, ProtOutputs):
            self.assertEqual(output, ProtOutput)

        # Sub-test 6: imprint ratio 1, sampling method 'sampled', ignore Cys and Trp
        # run carbonara 3 seq, inprint ratio 1, sampling method samlpled, ignored aminoacid Cys Trp
        args = {'atomStruct': PDB_output,
                'numSamples': 3,
                'imprintRatio' : 1,
                'bSampled' : 1,
                'ignoreAminoacids' : True,
                'selectIgnoredAminoacids' : "C, W"
                }
        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara inp_ratio 1 sampled C,W ignored')
        self.launchProtocol(protCarbonaraSamplingSequence)
        
        outPuts = []
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_0")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_1")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_2")

        ProtOutputs = []
        for t in protCarbonaraSamplingSequence.iterOutputAttributes():
            name = t[0]
            if name.startswith("DONOTSAVESESSION_Atom_struct__2_000042"):
                ProtOutputs.append(name)
        
        for output, ProtOutput in zip(outPuts, ProtOutputs):
            self.assertEqual(output, ProtOutput)

        # Sub-test 7: imprint ratio 1, sampling method 'sampled', ignore Cys/Trp, CPU mode
        # run carbonara 3 seq, inprint ratio 1, sampling method samlpled, ignored aminoacid Cys Trp, cpu used
        args = {'atomStruct': PDB_output,
                'numSamples': 3,
                'imprintRatio' : 1,
                'bSampled' : 1,
                'ignoreAminoacids' : True,
                'selectIgnoredAminoacids' : "C, W",
                'useGpu' : False
                }
        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara inp_ratio 1 sampled C,W ignored, cpu')
        self.launchProtocol(protCarbonaraSamplingSequence)
        
        outPuts = []
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_0")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_1")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_2")

        ProtOutputs = []
        for t in protCarbonaraSamplingSequence.iterOutputAttributes():
            name = t[0]
            if name.startswith("DONOTSAVESESSION_Atom_struct__2_000042"):
                ProtOutputs.append(name)
        
        for output, ProtOutput in zip(outPuts, ProtOutputs):
            self.assertEqual(output, ProtOutput)

        # Sub-test 8: imprint ratio 1, sampling method 'sampled', ignore Cys/Trp, CPU, AlphaFold
        # run carbonara 3 seq, inprint ratio 1, sampling method samlpled, 
        # ignored aminoacid Cys Trp, cpu used and launches AlphaFold
        #
        args = {'atomStruct': PDB_output,
                'numSamples': 3,
                'imprintRatio' : 1,
                'bSampled' : 1,
                'computeAlphaFold' : True,
                'ignoreAminoacids' : True,
                'selectIgnoredAminoacids' : "C, W",
                'useGpu' : False
                }
        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'carbonara inp_ratio 1 sampled C,W ignored, cpu')
        self.launchProtocol(protCarbonaraSamplingSequence)
        
        outPuts = []
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_0")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_1")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_2")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_0_binder")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_1_binder")
        outPuts.append(os.path.splitext(os.path.basename(PDB_output.getFileName()))[0] + "_2_binder")
        print(outPuts)
        ProtOutputs = []
        for t in protCarbonaraSamplingSequence.iterOutputAttributes():
            name = t[0]
            if name.startswith("DONOTSAVESESSION_Atom_struct__2_000042"):
                ProtOutputs.append(name)
        
        for output, ProtOutput in zip(outPuts, ProtOutputs):
            self.assertEqual(output, ProtOutput)
