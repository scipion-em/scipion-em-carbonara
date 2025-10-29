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
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from tkinter.messagebox import showerror
from pyworkflow.gui.text import openTextFileEditor

class CarbonaraViewer(ProtocolViewer):
    """ Visualize the output of protocol carbonara sampling sequence """
    _label= 'carbonara viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [CarbonaraSamplingSequence]

    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization of CARBonAra results')
        group = form.addGroup('Scores')
        group.addParam('displayScores', params.LabelParam,
                       label='CARBonAra sequence scores', 
                       help='List of ID sequences and CARBonAra scores')
        group = form.addGroup('Sequences')
        group.addParam('displaySequences', params.LabelParam,
                       label='CARBonAra sampled sequences', 
                       help='Sampled sequences in .fasta format. The starting sequence' \
                       ' is also included.')
        group = form.addGroup('Alignment')
        group.addParam('displayAlignment', params.LabelParam,
                       label='ClustalOmega alignment', 
                       help='Sampled sequences in .aln format. The starting sequence' \
                       ' is also included.')
        group = form.addGroup('Alignment Summary')
        group.addParam('displayAlignmentSummary', params.LabelParam,
                       label='Clustal Omega alignment summary', 
                       help='Starting sequence and ClustalOmega identity and similarity' \
                       ' symbols.')


    def _getVisualizeDict(self):
        return {
            'displayScores': self._showScores,
            'displaySequences': self._showSequences,
            'displayAlignment': self._showAlignment,
            'displayAlignmentSummary': self._showAlignmentSummary
        }
    
    # view all text files located in protocol extra folder

    def _showScores(self, obj, **args):
        self.dir_path = os.path.abspath(self.protocol._getExtraPath())
        # subdirectory = "carbonara_results"
        # self.subdir_path = os.path.join(self.dir_path, subdirectory)
        seqAndScoreFile = os.path.join(self.dir_path, "sorted_scores.txt")
        if os.path.exists(seqAndScoreFile ):
            openTextFileEditor(seqAndScoreFile)

    def _showSequences(self, obj, **args):
        fastaFile = os.path.join(self.dir_path, "merged.fasta")
        if os.path.exists(fastaFile):
            openTextFileEditor(fastaFile)

    def _showAlignment(self, obj, **args):      
        clustalAlignFile = os.path.join(self.dir_path, "clustal.aln")
        if os.path.exists(clustalAlignFile):
            openTextFileEditor(clustalAlignFile)

    def _showAlignmentSummary(self, obj, **args):  
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