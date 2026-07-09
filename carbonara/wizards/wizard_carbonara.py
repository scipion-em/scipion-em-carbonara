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
Wizard registrations for the CarbonaraSamplingSequence protocol.

Wizards provide interactive GUI helpers that let users select chains or
residues from a visualised structure and write the selection back into
the corresponding protocol parameter.
"""

from ..protocols import CarbonaraSamplingSequence
from pwem.wizards import SelectChainSWizard, SelectResidueWizard

# Chain-selection wizard: reads atomStruct, writes selectStructureChains
SelectChainSWizard().addTarget(protocol=CarbonaraSamplingSequence,
                               targets=['selectStructureChains'],
                               inputs=['atomStruct'],
                               outputs=['selectStructureChains'])

# Residue-selection wizard: currently disabled because the corresponding
# protocol parameters (selectKnownStructureResidues,
# selectUnknownStructureResidues) are commented out. Re-enable when restored.
# SelectResidueWizard().addTarget(protocol=CarbonaraSamplingSequence,
#                                 targets=['selectKnownStructureResidues',
#                                          'selectUnknownStructureResidues'],
#                                 inputs=['atomStruct'],
#                                 outputs=['selectKnownStructureResidues',
#                                          'selectUnknownStructureResidues'])
