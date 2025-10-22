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

from ..protocols import CarbonaraSamplingSequence
from pyworkflow.tests import BaseTest, setupTestProject
from pwem.protocols.protocol_import import (ProtImportPdb,
     ProtImportVolumes)


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)


class TestImportData(TestImportBase):
    """ Import atomic structure (PDBx/mmCIF file)
    """
    def _importStructurePDB(self, pdbID="5ni1"):  # Haemoglobine
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': pdbID
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import structure\n %s' % pdbID)
        self.launchProtocol(protImportPDB)
        self.assertTrue(protImportPDB.outputPdb.getFileName())
        return protImportPDB.outputPdb


class TestCarbonaraSamplingSequence(TestImportData):
    """ Test Carbonara protocol Sampling Sequence
    """

    def testCarbonara1(self):
        """ This test checks that carbonara runs with default parameters
            except number of sequences (3)
        """
        print("Run Carbonara Sampling Sequence protocol from imported pdb file " \
              "with default parameters and 3 predicted sequences")

        # import PDB
        Haemoglobin_PDB = self._importStructurePDB()
        self.assertTrue(Haemoglobin_PDB.getFileName())

        args = {'atomStruct': Haemoglobin_PDB,
                'numSamples': 3
                }

        protCarbonaraSamplingSequence = self.newProtocol(
            CarbonaraSamplingSequence, **args)
        protCarbonaraSamplingSequence.setObjLabel(
            'Haemoglobin 3 seqs default params')
        self.launchProtocol(protCarbonaraSamplingSequence)
        outPuts = ['5ni1_0', '5ni1_1', '5ni1_2']
        for output, ProtOutput in zip(outPuts, 
                                      protCarbonaraSamplingSequence.
                                      iterOutputAttributes()):
            self.assertTrue(str(ProtOutput), output)
