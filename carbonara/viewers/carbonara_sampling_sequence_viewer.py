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
import re
import pyworkflow.protocol.params as params
from ..protocols import CarbonaraSamplingSequence
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from tkinter.messagebox import showerror
from pyworkflow.gui.text import openTextFileEditor
from pwem.viewers import TableView, Chimera
from tkinter import messagebox
from pwem import Domain


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
        self.complex_file = None
        self.binder_file = None

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

            self.atom_structure_complex = []
            self.atom_structure_binder = []

            self.dir_path = os.path.abspath(self.protocol._getExtraPath())
            if os.path.exists(
                os.path.join(self.dir_path, "alphafold_predictions_multimer")):
                self.protocol.multimer_folder_path = \
                    os.path.join(self.dir_path, "alphafold_predictions_multimer")

                for filename in sorted(
                    os.listdir(self.protocol.multimer_folder_path)):
                    if filename.endswith(".pdb") or filename.endswith(".cif"):
                        self.atom_structure_complex.append(filename)

                group = form.addGroup('Complex structure viewer')
                group.addParam('displayComplexAtomStructure', params.EnumParam, 
                               display=params.EnumParam.DISPLAY_LIST,
                               choices=self.atom_structure_complex,
                               default = 0,
                               label='Select one COMPLEX', 
                               help='Select one complex with more than one chain '
                                    'to visualize its atom structure. If a binder ' \
                                    'was associated to this complex, the structure' \
                                    ' of the binder will be also shown.')
            else:

                if os.path.exists(
                    os.path.join(self.dir_path, "alphafold_predictions_binder")): 
                    self.protocol.binder_folder_path = \
                        os.path.join(self.dir_path, "alphafold_predictions_binder")           
            
                    for filename in sorted(
                        os.listdir(self.protocol.binder_folder_path)):
                        if filename.endswith(".pdb") or filename.endswith(".cif"):
                            self.atom_structure_binder.append(filename)

                    group = form.addGroup('Binder structure viewer')
                    group.addParam('displayBinderAtomStructure', params.EnumParam, 
                                   display=params.EnumParam.DISPLAY_LIST,
                                   choices=self.atom_structure_binder,
                                   default = 0,
                                   label='Select one BINDER',
                                   help='Select one binder (only one chain) to visualize'
                                        ' its atom structure')


    def _getVisualizeDict(self):
        vis_dict = {
            'displayScores': self._showAllScores,
            'displaySequences': self._showSequences,
            'displayAlignment': self._showAlignment,
            'displayAlignmentSummary': self._showAlignmentSummary,
        }

        # Add complex if it exits and binder does or does not
        if (self.protocol.computeAlphaFold and 
            os.path.exists(os.path.join(self.dir_path, "alphafold_predictions_multimer"))):
            
            self.multimer_folder_path = \
                os.path.join(self.dir_path, "alphafold_predictions_multimer")
            
            choices = self.atom_structure_complex
            idx = self.displayComplexAtomStructure.get()
            if idx is not None and 0 <= idx < len(choices):
                complex_choice = choices[idx]
                self.complex_file = os.path.join(self.multimer_folder_path,
                                    complex_choice)

            vis_dict['displayComplexAtomStructure'] = \
                lambda p=None: self._showStructuresInChimeraX([self.complex_file])

        # Add binder if it exits and complex does not
        if (self.protocol.computeAlphaFold and not
            os.path.exists(os.path.join(self.dir_path, "alphafold_predictions_multimer")) and
            os.path.exists(os.path.join(self.dir_path, "alphafold_predictions_binder"))):

            self.binder_folder_path = \
                os.path.join(self.dir_path, "alphafold_predictions_binder")
            
            choices = self.atom_structure_binder
            idx = self.displayBinderAtomStructure.get()
            
            if idx is not None and 0 <= idx < len(choices):
                binder_choice = choices[idx]
                self.binder_file = os.path.join(self.binder_folder_path,
                                   binder_choice )

            vis_dict['displayBinderAtomStructure'] = \
                lambda p=None: self._showStructuresInChimeraX([self.binder_file])
        """
        TODO: Comment
        # Add complex and binder if they both exist
        if (self.protocol.computeAlphaFold and 
            os.path.exists(os.path.join(self.dir_path, "alphafold_predictions_multimer")) and 
            os.path.exists(os.path.join(self.dir_path, "alphafold_predictions_binder"))):
            
            self.multimer_folder_path = \
                os.path.join(self.dir_path, "alphafold_predictions_multimer")
            self.binder_folder_path = \
                os.path.join(self.dir_path, "alphafold_predictions_binder")
            
            # Complex
            choices = self.atom_structure_complex
            idx = self.displayComplexAtomStructure.get()
            complex_file = None
            if idx is not None and 0 <= idx < len(choices):
                complex_choice = choices[idx]
                complex_file = os.path.join(self.multimer_folder_path,
                                    complex_choice)
                self.complex_file = complex_file
                print("self.complex_file: ", self.complex_file)

            
            # Binder
            choices = self.atom_structure_binder
            idx = self.displayBinderAtomStructure.get()
            binder_file = None
            if idx is not None and 0 <= idx < len(choices):
                binder_choice = choices[idx]
                self.binder_file = os.path.join(self.binder_folder_path,
                                   binder_choice )
                self.binder_file = binder_file
                print("self.binder_file: ", self.binder_file)
            
            if complex_file and binder_file:
                vis_dict['displayComplexAndBinder'] = \
                lambda p=None: self._showStructuresInChimeraX([self.complex_file, self.binder_file])
        """    
        return vis_dict


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
        self.dir_path = os.path.abspath(self.protocol._getExtraPath())
        csvFile = os.path.join(self.dir_path, "sorted_scores.csv")

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

    def _showStructuresInChimeraX(self, e=None):

        # Atom structures to show
        fileNames = []
        if self.complex_file:
            fileNames.append(self.complex_file)
            if (os.path.exists(os.path.join(self.dir_path, "alphafold_predictions_binder"))):
                self.binder_folder_path = \
                    os.path.join(self.dir_path, "alphafold_predictions_binder")
                baseFileName = self.complex_file.split("/")[-1].split("_unrelaxed_")[0]

                for filename in sorted(
                        os.listdir(self.binder_folder_path)):
                        if (filename.startswith(baseFileName + "_chain") and
                            filename.endswith(".pdb" or ".cif")):
                            binder_file = os.path.join(self.binder_folder_path, filename)
                            fileNames.append(binder_file)
                
        if self.binder_file:
            fileNames.append(self.binder_file)

        fileNames = [os.path.abspath(f) for f in fileNames if f]

        bildFileName = self.protocol._getExtraPath("axis_output.bild")
        fnCmd = self.protocol._getExtraPath("chimera_output.cxc")
        
        # Create axes
        dim = 150.
        sampling = 1.

        Chimera.createCoordinateAxisFile(
            dim,
            bildFileName=bildFileName,
            sampling=sampling
        )
         # Write ChimeraX script
        
        with open(fnCmd, 'w', encoding='utf-8') as f:
            # change to workingDir
            # If we do not use cd and the project name has an space
            # the protocol fails even if we pass absolute paths
            f.write('cd %s\n' % os.getcwd())
            # reference axis model = 1
            f.write("open %s\n" % bildFileName)
            f.write("cofr 0,0,0\n")  # set center of coordinates
            # open all fileNames
            for pdbFileName in fileNames:
                f.write("open %s\n" % pdbFileName)
            f.write("color by bfactor palette alphafold\n")
            if len(fileNames) > 1:
                try:
                    filename = fileNames[1]
                    filename = os.path.basename(pdbFileName)
                    match = re.search(r'chain([A-Za-z0-9]+)', filename)
                    if match:
                        chain = match.group(1)
                        f.write("mmaker #3/%s to #2/%s\n" % (chain, chain))
                except Exception as e:
                    print("Error writing mmaker command:", e)

        # run in the background
        chimeraPlugin = Domain.importFromPlugin('chimera', 'Plugin', doRaise=True)
        chimeraPlugin.runChimeraProgram(chimeraPlugin.getProgram(), fnCmd + "&",
                                        cwd=os.getcwd())
        return []

    def errorWindow(tkParent, msg):
        try:
            showerror("Error",  # bar title
                  msg,  # message
                  parent=tkParent)
        except Exception as ex:
            print("Error:", msg)
