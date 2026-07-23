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

"""
Viewer for the CarbonaraSamplingSequence protocol.

Provides GUI access to:
  - CARBonAra sequence scores (table view)
  - Sampled sequences (FASTA text editor)
  - ClustalO alignment (text editor)
  - ClustalO alignment summary (text editor)
  - AlphaFold-predicted structures (ChimeraX viewer, if computed)
"""

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
    """Show an error dialog, falling back to print if Tk is unavailable.

    The tkParent parameter ensures the dialog appears above other windows.
    """
    try:
        # if tkRoot is null the error message may be behind
        # other windows
        messagebox.showerror("Error",  # bar title
                              msg,  # message
                              parent=tkParent)
    except:
        print(("Error:", msg))


class CarbonaraViewer(ProtocolViewer):
    """Visualize the output of the CarbonaraSamplingSequence protocol.

    Provides four always-available sections:
      - Scores: table of CARBonAra and AlphaFold metrics
      - Sequences: merged FASTA file with all sampled sequences
      - Alignment: full ClustalO alignment
      - Alignment Summary: filtered alignment showing reference + consensus

    If AlphaFold was computed, additional sections appear:
      - Complex structure viewer (multi-chain predictions)
      - Binder structure viewer (single-chain predictions)
    Structures are opened in ChimeraX with alphafold B-factor coloring.
    """
    _label= 'carbonara viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [CarbonaraSamplingSequence]

    def __init__(self, **kwargs):
        """Initialise the viewer and reset structure file trackers."""
        ProtocolViewer.__init__(self, **kwargs)
        self.complex_file = None
        self.binder_file = None

    def _defineParams(self, form):
        """Define the viewer form parameters.

        Always shows: Scores, Sequences, Alignment, Alignment Summary.
        Conditionally shows: Complex structure viewer or Binder structure viewer
        (depending on which AlphaFold outputs exist).
        """
        form.addSection(label='Visualization of CARBonAra results')

        # --- Always-available sections ------------------------------------
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
        
        self.dir_path = os.path.abspath(self.protocol._getExtraPath())
        
        # --- Conditional AlphaFold structure sections ---------------------
        if self.protocol.computeAlphaFold:

            self.atom_structure_complex = []
            self.atom_structure_binder = []

            multimer_folder_path = os.path.join(self.dir_path, "alphafold_predictions_multimer")
            binder_folder_path   = os.path.join(self.dir_path, "alphafold_predictions_binder")

            multimer_exists = os.path.exists(multimer_folder_path)
            binder_exists   = os.path.exists(binder_folder_path)

            # Show complex viewer if multimer predictions exist
            if multimer_exists:
                for filename in sorted(
                    os.listdir(multimer_folder_path)):
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
            
            # Show binder viewer if binder predictions exist
            if binder_exists:
                for filename in sorted(
                    os.listdir(binder_folder_path)):
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

            # Show complex and binder viewer if both exist
            if multimer_exists and binder_exists:
                group = form.addGroup('Complex + Binder structure viewer')
                group.addParam('displayComplexAndBinder', params.BooleanParam, 
                                default = False,
                                label='Show COMPLEX and BINDER together',
                                help='Visualize both structures simultaneously in ChimeraX')

    def _getVisualizeDict(self):
        """Return a mapping from parameter names to visualization callbacks.

        Scipion calls this method to determine which function to invoke when
        the user clicks a "Show" button in the viewer GUI. Each callback
        opens the appropriate file or launches ChimeraX for structures.
        """
        vis_dict = {
            'displayScores': self._showAllScores,
            'displaySequences': self._showSequences,
            'displayAlignment': self._showAlignment,
            'displayAlignmentSummary': self._showAlignmentSummary,
        }

        # --- Detect AlphaFold outputs once ---------------------------------

        multimer_folder_path = os.path.join(self.dir_path, "alphafold_predictions_multimer")
        binder_folder_path   = os.path.join(self.dir_path, "alphafold_predictions_binder")

        multimer_exists = os.path.exists(multimer_folder_path)
        binder_exists   = os.path.exists(binder_folder_path)

        # Add complex structure viewer if multimer predictions exist
        if multimer_exists:
            self.multimer_folder_path = multimer_folder_path
            
            idx = self.displayComplexAtomStructure.get()
            if idx is not None and 0 <= idx < len(self.atom_structure_complex):
                complex_choice = self.atom_structure_complex[idx]
                self.complex_file = os.path.join(self.multimer_folder_path,
                                    complex_choice)

            vis_dict['displayComplexAtomStructure'] = \
                lambda p=None: self._showStructuresInChimeraX([self.complex_file])

        # Add binder structure viewer if binder predictions exist 
        if binder_exists:
            self.binder_folder_path = binder_folder_path
            
            idx = self.displayBinderAtomStructure.get()            
            if idx is not None and 0 <= idx < len(self.atom_structure_binder):
                binder_choice = self.atom_structure_binder[idx]
                self.binder_file = os.path.join(self.binder_folder_path,
                                   binder_choice )

            vis_dict['displayBinderAtomStructure'] = \
                lambda p=None: self._showStructuresInChimeraX([self.binder_file])
            
        # Add complex and binder structure viewer if both exist 
        if multimer_exists and binder_exists:
            self.binder_folder_path = binder_folder_path

            vis_dict['displayComplexAndBinder'] = \
                lambda p=None: self._showStructuresInChimeraX(
                    [self.complex_file, self.binder_file]
                )    
        
        return vis_dict


    def _showScores(self, headerList, dataList, mesg, title):
        """Display a table of scores using Scipion's TableView widget.

        Shows an error dialog if dataList is empty.
        """
        if not dataList:
            errorWindow(self.getTkRoot(), "No data available")
            return

        TableView(headerList=headerList,
                  dataList=dataList,
                  mesg=mesg,
                  title=title,
                  height=len(dataList), width=250, padding=40)

    def _showAllScores(self, e=None):
        """Load the sorted_scores.csv file and display it in a table.

        The 'sequence' column is removed from the display to keep the table
        compact (sequences can be viewed in the FASTA file instead).
        """
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
        """Open the merged FASTA file in a text editor."""
        fastaFile = os.path.join(self.dir_path, "merged.fasta")
        if os.path.exists(fastaFile):
            openTextFileEditor(fastaFile)

    def _showAlignment(self, obj, **args):
        """Open the full ClustalO alignment file in a text editor."""
        clustalAlignFile = os.path.join(self.dir_path, "clustal.aln")
        if os.path.exists(clustalAlignFile):
            openTextFileEditor(clustalAlignFile)

    def _showAlignmentSummary(self, obj, **args):
        """Open the filtered ClustalO alignment summary in a text editor.

        The summary shows only the reference sequence and consensus symbols.
        """
        clustalAlignFileSummaryFile = os.path.join(self.dir_path, "clustal_summary.aln")
        if os.path.exists(clustalAlignFileSummaryFile):
            openTextFileEditor(clustalAlignFileSummaryFile)

    def _showStructuresInChimeraX(self, fileNames, e=None):
        """Launch ChimeraX to visualize AlphaFold-predicted structures.

        Generates a ChimeraX script (.cxc) that:
          - Opens coordinate axes for reference
          - Opens all selected structure files (complex and/or binder)
          - Colors by B-factor using the AlphaFold pLDDT palette
          - Optionally superposes binder onto complex using MatchMaker

        The script is written to chimera_output.cxc and executed in the
        background via the chimera plugin.
        """
        # Normalize and filter
        fileNames = [os.path.abspath(f) for f in fileNames if f]

        if not fileNames:
            print("No structures to visualize.")
            return []
        
        # Output script paths
        bildFileName = self.protocol._getExtraPath("axis_output.bild")
        fnCmd = self.protocol._getExtraPath("chimera_output.cxc")

        # Create coordinate axes BILD file
        Chimera.createCoordinateAxisFile(
            dim=150.,
            bildFileName=bildFileName,
            sampling=1.
        )
        
         # Write ChimeraX script
        
        with open(fnCmd, 'w', encoding='utf-8') as f:
            # change to workingDir (avoid issues with spaces)
            f.write('cd %s\n' % os.getcwd())

            # reference axis model = 1
            f.write("open %s\n" % bildFileName)
            f.write("cofr 0,0,0\n")  # set center of coordinates

            # open all structure files
            for pdbFileName in fileNames:
                f.write("open %s\n" % pdbFileName)

            # Color by B-factor using AlphaFold pLDDT palette
            f.write("color by bfactor palette alphafold\n")

            # If multiple structures, superpose binder onto complex using MatchMaker
            if len(fileNames) > 1:
                try:
                    filename = os.path.basename(fileNames[1])
                    # Extract chain ID from filename (e.g. "seq1_chainA.pdb" -> "A")
                    match = re.search(r'chain([A-Za-z0-9]+)', filename)
                    if match:
                        chain = match.group(1)
                        print("chain: ", chain)
                        # Superpose chain X of model #3 onto chain X of model #2
                        f.write("mmaker #3/%s to #2/%s\n" % (chain, chain))
                except Exception as e:
                    print("Error writing mmaker command:", e)

        # Launch ChimeraX in the background with the generated script
        chimeraPlugin = Domain.importFromPlugin('chimera', 'Plugin', doRaise=True)
        chimeraPlugin.runChimeraProgram(chimeraPlugin.getProgram(), fnCmd + "&",
                                        cwd=os.getcwd())
        return []
    