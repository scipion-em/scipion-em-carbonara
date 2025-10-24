# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
#                Roberto Marabini (roberto@cnb.csic.es)
# *
# * CNB CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
import pyworkflow.protocol.params as params
from ..protocols import CarbonaraSamplingSequence
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer, WEB_DJANGO
from tkinter.messagebox import showerror
from pyworkflow.gui.text import openTextFileEditor

SEQANDSCORE_FILE = 0
FASTA_FILE= 1
CLUSTAL_FILE = 2
CLUSTALSUMM_FILE = 3

class CarbonaraViewer(Viewer):
    """ Visualize the output of protocol carbonara sampling sequence """
    _label= 'carbonara viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [CarbonaraSamplingSequence]

    def __init__(self, **kwargs):
        Viewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization of CARBonAra result summary text files')
        form.addParam('displayFiles', params.EnumParam,
                      choices=['Sequences with scores file', 'Sequences in fasta format', 
                               "Alignment with Clustal Omega", "Summary of the Clustal Alignment"],
                      default=SEQANDSCORE_FILE,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Select the file to visualize',
                      help='')

    def _getVisualizeDict(self):
        return {
            'displayFiles': self._showTextFiles,
        }
    
    def _showTextFiles(self, obj, **args):
        self.dir_path = os.path.abspath(self.protocol._getExtraPath())
        subdirectory = "carbonara_results"
        self.subdir_path = os.path.join(self.dir_path, subdirectory)
        
        # view all text files located in protocol extra folder

        if self.displayFiles == SEQANDSCORE_FILE:
            seqAndScoreFile = os.path.join(self.dir_path, "sorted_scores.txt")
            if os.path.exists(seqAndScoreFile ):
                openTextFileEditor(seqAndScoreFile)

        if self.displayFiles == FASTA_FILE:
            fastaFile = os.path.join(self.dir_path, "merged.fasta")
            if os.path.exists(fastaFile):
                openTextFileEditor(fastaFile)

        if self.displayFiles == CLUSTAL_FILE:
            clustalAlignFile = os.path.join(self.dir_path, "clustal.aln")
            if os.path.exists(clustalAlignFile):
                openTextFileEditor(clustalAlignFile)

        if self.displayFiles == CLUSTALSUMM_FILE:
            clustalAlignFileSummaryFile = os.path.join(self.dir_path, "clustal_summary.aln")
            if os.path.exists(clustalAlignFileSummaryFile):
                openTextFileEditor(clustalAlignFileSummaryFile)
    
    def errorWindow(tkParent, msg):
        try:
            showerror("Error",  # bar title
                  msg,  # message
                  parent=tkParent)
        except Exception as ex:
            print("Error:", msg)