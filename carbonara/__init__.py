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
Plugin entry point for scipion-em-carbonara.

Defines the Plugin class which:
  - registers the CARBONARA_ENV_ACTIVATION variable,
  - builds the conda-activation command used to invoke the carbonara binary,
  - installs the required conda environments (carbonara + colabfold) via
    defineBinaries.
"""
import os
import sys

import pwem
import pyworkflow.utils as pwutils
from scipion.install.funcs import VOID_TGZ
import json
# from pwem.constants import MAXIT
import subprocess

from .constants import (CARBONARA_ENV_ACTIVATION, conda_env,
                        GIT_CLONE_CMD, CARBONARA_PROGRAM,
                        colab_env, colabfold_repo, jax_api)

__version__ = "0.0.5"
_logo = "icon.png"
_references = ['Krapp2024']


class Plugin(pwem.Plugin):
    """Scipion plugin wrapper for the CARBonAra protein-design tool.

    Manages two conda environments:
      * carbonara  - runs the CARBonAra CLI and clustalo
      * colabfold  - runs ColabFold / AlphaFold for optional scoring
    """

    @classmethod
    def _defineVariables(cls):
        """Register the CARBONARA_ENV_ACTIVATION variable in Scipion's
        configuration so users can override the activation command."""
        cls._defineVar(CARBONARA_ENV_ACTIVATION,
                       cls.getCARBonAraActivationCmd())

    @classmethod
    def getCarbonaraCmd(cls):
        """Return the full shell command to invoke the carbonara binary.

        The command chains conda activation with the carbonara executable
        so that the correct environment is active at runtime.
        """
        cmd = cls.getCondaActivationCmd()
        cmd += " "
        cmd += cls.getVar(CARBONARA_ENV_ACTIVATION)
        cmd += f" && {CARBONARA_PROGRAM} "
        return cmd

    @classmethod
    def getEnviron(cls):
        """Build an Environ with the carbonara conda-env bin/ on PATH.

        This is needed so that clustalo (installed inside the carbonara
        conda env) can be found by alignment routines.
        Queries 'conda env list --json' to locate the env directory.
        """
        environ = pwutils.Environ(os.environ)
        # Locate the carbonara conda env and prepend its bin/ to PATH
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
    def getCARBonAraActivationCmd(cls):
        """Return the conda activation command for the carbonara env."""
        return f'conda activate {conda_env}'

    @classmethod
    def getDependencies(cls):
        """Return a list of dependencies. Include conda if
        activation command was not found."""
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

    # @classmethod
    def defineBinaries(self, env):
        """Register the carbonara package with Scipion's binary manager.

        Two conda envs are created/updated:
          1. carbonara  - python 3.9, clustalo, CARBonAra from pinned commit
          2. colabfold  - python 3.10, jax with CUDA, ColabFold with alphafold
        A flag file marks successful installation so Scipion can skip
        redundant reinstalls.
        """
        def create_or_update_conda_env(
                env_name, packages=None, python_version="3.8"):
            """Return a conda shell command that creates env_name if it
            does not exist, or updates it otherwise.

            If packages is None and the env already exists, a full
            'conda update --all' is issued instead.
            """
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
                # command = f"conda install --yes --name {env_name}"
                # if packages:
                #     command += packages
                # else:
                command = f"conda update --all --yes --name {env_name}"
            else:
                print(f"Creating new Conda environment: {env_name}")
                command = "conda create --yes --name "\
                          f"{env_name} python={python_version}"
                # if packages:
                #     command += packages

            return command

        def getCARBonAraInstallation(version):
            """Build the install command and register it with Scipion.

            The command chains env creation, clustalo install, CARBonAra
            pip-install, colabfold env creation, jax/colabfold install,
            and a final 'touch' to write the flag file.
            """
            FLAG = f'{conda_env}_{version}_installed'

            conda_env_path = os.environ.get('CONDA_PREFIX', sys.prefix)
            conda_env_path = os.path.dirname(conda_env_path)
            conda_env_path = os.path.join(conda_env_path, conda_env)

            FLAG = os.path.join(conda_env_path, FLAG)

            install_cmd = [
                # Activate conda and create the carbonara environment
                f'{cls.getCondaActivationCmd()}'
                f'{create_or_update_conda_env(conda_env, python_version="3.9")} && '
                # Install clustalo inside the carbonara env
                f'conda run -n {conda_env} conda install -y -c bioconda clustalo && '
                # Install CARBonAra from the pinned git commit
                f'conda run -n {conda_env} {GIT_CLONE_CMD} && '
                # Create the colabfold environment
                f'{create_or_update_conda_env(colab_env, python_version="3.10")} && '
                # Install jax with CUDA support
                f'conda run -n {colab_env} pip install --upgrade "jax[cuda12]" -f {jax_api} && '
                # Install ColabFold with alphafold extras
                f'conda run -n {colab_env} pip install "{colab_env}[alphafold] @ {colabfold_repo}" && '
                # Write the flag file so Scipion knows installation succeeded
                f'touch {FLAG}'
            ]

            finalCmds = [(install_cmd, FLAG)]
            print("finalCmds: ", finalCmds)

            # Register the package with Scipion's binary manager
            env.addPackage(
                name="carbonara",
                version=version,
                tar=VOID_TGZ,  # VOID_TGZ is a no-op tar; the real work is in commands
                commands=finalCmds,
                neededProgs=cls.getDependencies(),
                default=True)

        getCARBonAraInstallation(version=__version__)

