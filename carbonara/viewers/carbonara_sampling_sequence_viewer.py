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
import csv
import pyworkflow.protocol.params as params
from ..protocols import CarbonaraSamplingSequence
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from tkinter.messagebox import showerror
from pyworkflow.gui.text import openTextFileEditor
from pwem.viewers import TableView
from tkinter import messagebox

def errorWindow(tkParent, msg):
    try:
        # if tkRoot is null the error message may be behind
        # other windows
        messagebox.showerror("Error",  # bar title
                              msg,  # message
                              parent=tkParent)
    except:
        print(("Error:", msg))

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
        
        if self.protocol.computeAlphaFold:
            atom_structure_complex = []
            atom_structure_binder = []
            dir_path = os.path.abspath(self.protocol._getExtraPath())

            if self.protocol.COMPLEX == True:
                for filename in sorted(
                    os.listdir(self.protocol.multimer_folder_path)):
                    if filename.endswith(".pdb") or filename.endswith(".cif"):
                        atom_structure_complex.append(filename)

                group = form.addGroup('Complex structure viewer')
                group.addParam('displayAtomStructure', 
                       choices=atom_structure_complex,
                       display=params.EnumParam.DISPLAY_LIST,
                       label='Select one complex to visualize its atom structure', 
                       help='')
                        
            if self.BINDER == True and len(self.chains) == 1:
                for filename in sorted(
                    os.listdir(self.protocol.binder_folder_path)):
                    if filename.endswith(".pdb") or filename.endswith(".cif"):
                        atom_structure_binder.append(filename)

                group = form.addGroup('Binder structure viewer')
                group.addParam('displayAtomStructure', 
                       choices=atom_structure_binder,
                       display=params.EnumParam.DISPLAY_LIST,
                       label='Select one binder to visualize its atom structure', 
                       help='')
                       
        

    def _getVisualizeDict(self):
        return {
            'displayScores': self._showAllScores,
            'displaySequences': self._showSequences,
            'displayAlignment': self._showAlignment,
            'displayAlignmentSummary': self._showAlignmentSummary
        }
    
    # view all text files located in protocol extra folder

    def _showScores(self, headerList, dataList, mesg, title):

        if not dataList:
            errorWindow(self.getTkRoot(), "No data available")
            return

        TableView(headerList=headerList,
                  dataList=dataList,
                  mesg=mesg,
                  title=title,
                  height=len(dataList), width=250, padding=40)

    def _showAllScores(self, e=None):
        dir_path = os.path.abspath(self.protocol._getExtraPath())
        csvFile = os.path.join(dir_path, "sorted_scores.csv")

        headers = []
        dataList = []
        
        with open(csvFile, newline="", encoding='utf-8') as csvF:
            csvReader = csv.reader(csvF)
            headers = next(csvReader)
        
            # remove sequence column to simplify the table output
            if "sequence" in headers:
                idx = headers.index("sequence")
                headers.pop(idx)

                dataList = []
                for row in csvReader:
                    row.pop(idx)
                    dataList.append(row)
            else:
                dataList = list(csvReader)

        self.headerList = headers
        self.dataList = dataList  

        mesg = ""
        title = "Table of CARBonAra and AlphaFold scores"
        self._showScores(headers, dataList, mesg, title)

    def _showSequences(self, obj, **args):
        self.dir_path = os.path.abspath(self.protocol._getExtraPath())
        fastaFile = os.path.join(self.dir_path, "merged.fasta")
        if os.path.exists(fastaFile):
            openTextFileEditor(fastaFile)

    def _showAlignment(self, obj, **args): 
        self.dir_path = os.path.abspath(self.protocol._getExtraPath())     
        clustalAlignFile = os.path.join(self.dir_path, "clustal.aln")
        if os.path.exists(clustalAlignFile):
            openTextFileEditor(clustalAlignFile)

    def _showAlignmentSummary(self, obj, **args): 
        self.dir_path = os.path.abspath(self.protocol._getExtraPath())
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