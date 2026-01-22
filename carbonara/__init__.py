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
import sys

import pwem
import pyworkflow.utils as pwutils
from scipion.install.funcs import VOID_TGZ
import subprocess
import json
# from pwem.constants import MAXIT

from .constants import (CARBONARA_ENV_ACTIVATION, conda_env,
                        GIT_CLONE_CMD, CARBONARA_PROGRAM,
                        colab_env, colabfold_repo, jax_api)

__version__ = "0.0.5"
_logo = "icon.png"
_references = ['Krapp2024']


class Plugin(pwem.Plugin):

    @classmethod
    def _defineVariables(cls):
        # CARBONARA_ENV_ACTIVATION="conda activate carbonara"
        cls._defineVar(CARBONARA_ENV_ACTIVATION,
                       cls.getCARBonAraActivationCmd())

    @classmethod
    def getCarbonaraCmd(cls):
        cmd = cls.getCondaActivationCmd()
        cmd += " "
        cmd += cls.getVar(CARBONARA_ENV_ACTIVATION)
        cmd += f" && {CARBONARA_PROGRAM} "
        return cmd

    @classmethod
    def getEnviron(cls):
        environ = pwutils.Environ(os.environ)
        # add conda bin path since clustalo has been installed here
        result = subprocess.run(
            ["conda", "env", "list", "--json"],
             capture_output=True,
             text=True,
             check=True
        )
        envs = json.loads(result.stdout)["envs"]
        env_path = next(p for p in envs if p.endswith(f"/{conda_env}"))
        bin_dir = env_path + ("/Scripts" if os.name == "nt" else "/bin")

        environ.update({'PATH': bin_dir},
                       position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def getCARBonAraActivationCmd(cls):   # getActivationCmd
        return f'conda activate {conda_env}'

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
        import subprocess
        import json

        def create_or_update_conda_env(
                env_name, packages=None, python_version="3.8"):
            # Get list of existing environments
            result = subprocess.run(
                ["conda", "env", "list", "--json"],
                capture_output=True,
                text=True,
                check=True
            )
            envs = json.loads(result.stdout)["envs"]
            env_exists = any(
                env.endswith(f"/{env_name}") or
                env.endswith(f"\\{env_name}") for env in envs)

            if env_exists:
                print(f"Updating existing Conda environment: {env_name}")
                command = f"conda install --yes --name {env_name}"
                if packages:
                    command += packages
                else:
                    command = f"conda update --all --yes --name {env_name}"
            else:
                print(f"Creating new Conda environment: {env_name}")
                command = "conda create --yes --name "\
                          f"{env_name} python={python_version}"
                if packages:
                    command += packages

            return command

        def getCARBonAraInstallation(version):
            FLAG = f'{conda_env}_{version}_installed'

            conda_env_path = os.environ.get('CONDA_PREFIX', sys.prefix)
            conda_env_path = os.path.dirname(conda_env_path)
            conda_env_path = os.path.join(conda_env_path, conda_env)

            FLAG = os.path.join(conda_env_path, FLAG)

            install_cmd = [
                # Activate Conda y create environment for carbonara
                f'{cls.getCondaActivationCmd()}'
                f'{create_or_update_conda_env(conda_env, python_version="3.9")} && '
                # Install
                f'conda run -n {conda_env} conda install -y -c bioconda clustalo && '
                f'conda run -n {conda_env} {GIT_CLONE_CMD} && '
                f'{create_or_update_conda_env(colab_env, python_version="3.10")} && '
                f'conda run -n {colab_env} pip install --upgrade "jax[cuda12]" -f {jax_api} && '
                f'conda run -n {colab_env} pip install "{colab_env}[alphafold] @ {colabfold_repo}" && '
                f'touch {FLAG}'
            ]

            finalCmds = [(install_cmd, FLAG)]
            print("finalCmds: ", finalCmds)

            # CARBonAra package registered in Scipion environment
            env.addPackage(
                name="carbonara",
                version=version,
                tar=VOID_TGZ,  # Only marker
                commands=finalCmds,
                neededProgs=cls.getDependencies(),
                default=True)

        getCARBonAraInstallation(version=__version__)

