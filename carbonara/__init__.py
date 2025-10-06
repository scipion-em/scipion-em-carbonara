# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
#                Roberto Marabini (roberto@cnb.csic.es)
# *
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
import pyworkflow.utils as pwutils
import pwem
from pyworkflow import Config

from .constants import *

__version__ = '0.0.1'
_logo = "logo.png"
_references = ['Krapp2024']


class Plugin(pwem.Plugin):
    _url = "https://github.com/scipion-em/scipion-em-carbonara"
    _supportedVersions = v1  # binary version

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(CARBONARA_ENV_ACTIVATION)
        cls._defineEmVar(DEFAULT_ACTIVATION_CMD)   
        
    @classmethod
    def getCARBonAraEnvActivation(cls):
        """ Remove the scipion home and activate the conda environment. """
        activation = cls.getVar(CARBONARA_ENV_ACTIVATION)
        scipionHome = Config.SCIPION_HOME + os.path.sep

        return activation.replace(scipionHome, "", 1)      

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch CARBonAra. """
        environ = pwutils.Environ(os.environ)

        # ...

        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. Include conda if
        activation command was not found. """
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')
            
        return neededProgs

    @classmethod
    def defineBinaries(cls, env):
            cls.addCARBonAraPackage(env, version=__version__)
            
    @classmethod
    def addCARBonAraPackage(env, version=__version__):
    	# try to get CONDA activation command
        installCmds = [
            cls.getCondaActivationCmd(),
            f'conda create -n carbonara',
            f'conda activate carbonara',
            f'pip install -e . &&'
        ]
    
        url = "https://github.com/LBM-EPFL/CARBonAra"
        gitCmds = [
            'cd .. &&',
            f'git clone {url} &&',
            f'cd CARBonAra;'
        ]
        gitCmds.extend(installCmds)
        carbonaraCmds = [(" ".join(gitCmds))]
        env.addPackage('carbonara', version=version,
                       tar='void.tgz',
                       commands=carbonaraCmds,
                       neededProgs=cls.getDependencies(),
                       default=True)
                       
    @classmethod
    def getActivationCmd(cls):
        """ Return the activation command. """
        return '%s %s' % (cls.getCondaActivationCmd(),
                          cls.getCARBonAraEnvActivation())

    @classmethod
    def getProgram(cls, program, gpus='0'):
        """ Create CARBonAracommand line. """
        fullProgram = '%s && CUDA_VISIBLE_DEVICES=%s carbonara %s' % (
            cls.getActivationCmd(), gpus, program)

        return fullProgram
                       
                       
