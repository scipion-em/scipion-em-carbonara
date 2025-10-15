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
# from pyworkflow import Config
from scipion.install.funcs import VOID_TGZ

from .constants import (CARBONARA_HOME, CARBONARA_ENV_ACTIVATION, V1)

__version__ = '0.0.1'  # plugin version
_logo = "icon.png"
_references = ['Krapp2024']


class Plugin(pwem.Plugin):
    _homeVar = CARBONARA_HOME
    _pathVars = [CARBONARA_HOME, CARBONARA_ENV_ACTIVATION]
    _url = "https://github.com/scipion-em/scipion-em-carbonara"
    _supportedVersions = [V1]  # binary version

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(CARBONARA_HOME, 'carbonara')
        cls._defineEmVar(CARBONARA_ENV_ACTIVATION,
                         cls.getCARBonAraActivationCmd())

    @classmethod
    def getCARBonAraActivationCmd(cls):
        return "conda activate carbonara && "

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

    # @classmethod
    def defineBinaries(cls, env):
        """
        How to install CARBonAra package from github in Scipion environment.
        """
        def getCARBonAraInstallation(version):
            repo_url = "https://github.com/LBM-EPFL/CARBonAra"
            repo_name = "CARBonAra"
            conda_env = "carbonara"
            FLAG = f'{conda_env}_{version}_installed'

            # Comands to clone repo if the repo folder doesn't
            # exist and configure environment
            installation_path = os.path.join(
                pwem.Config.EM_ROOT, f"{conda_env}-{version}")
            if os.path.isdir(installation_path):
                command = f'git -C {repo_name} pull '
            else:
                # Cloning repo
                command = f'git clone {repo_url} {repo_name} '

            install_cmds = [
                f'{command} && '
                f'cd {repo_name}  && '
                # Activate Conda y create environment
                f'{cls.getCondaActivationCmd()}'
                f' conda create -y -n {conda_env} python=3.9 && '
                f'{cls.getCARBonAraActivationCmd()}'
                # Install
                f'pip install . && '
                # Flag installation finished
                f'cd .. && '
                f'touch {FLAG}'
            ]

            finalCmds = [(" ".join(install_cmds), FLAG)]
            print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA: ",
                  os.getcwd(), os.path.join(pwem.Config.EM_ROOT, f"{conda_env}-{version}"),
                  finalCmds)

            # CARBonAra package registered in Scipion environment
            env.addPackage(
                name="carbonara",
                version=version,
                tar=VOID_TGZ,  # Only marker
                commands=finalCmds,
                neededProgs=cls.getDependencies(),
                default=True)

        getCARBonAraInstallation(version=V1)

    @classmethod
    def getProgram(cls, program, gpus='0'):
        """ Create CARBonAracommand line. """
        fullProgram = '%s && CUDA_VISIBLE_DEVICES=%s carbonara %s' % (
            cls.getActivationCmd(), gpus, program)

        return fullProgram
