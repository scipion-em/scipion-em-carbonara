=======================
Carbonara scipion plugin
=======================

This plugin allows to use run **CARBonAra** commands within the **Scipion** framework.

CARBonAra is a deep learning framework that facilitates protein sequence design by leveraging atomic coordinates, allowing for context-aware sequence generation. This method is particularly useful for integrating protein design with molecular environments, including non-protein entities, providing more control to protein engineering (See `https://github.com/LBM-EPFL/CARBonAra/>`_  for details)

Features:

Geometric Transformer: The framework uses a geometric transformer model based only on atomic coordinates and atomic elements, allowing it to handle any protein backbone scaffolds and various molecular environments.
    
Context Awareness: CARBonAra's design accounts for molecular environments, including non-protein entities, providing context-aware sequence generation.
    
Imprint Sequence Sampling: CARBonAra's imprint sampling method provides diverse sequences, balancing design flexibility with high-confidence predictions.


==========================
Install this plugin
==========================

You will need to use `3.0.0 <https://scipion-em.github.io/docs/release-3.0.0/>`_ version of Scipion to run these protocols.
Instructions to install this plugin (2 options):

- **Stable version**  

.. code-block:: 

      scipion3 installp -p scipion-em-carbonara
      
      cd scipion-em-carbonara
      rm -rf .git
      git init
    
      
OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version**

1. Download repository:

.. code-block::

    git clone git@github.com:scipion-em/scipion-em-carbonara.git
2. Install:

.. code-block::

    cd scipion-em-carbonara
    
    scipion3 installp -p path_to_scipion-em-carbonara --devel
    
- **Binary files**

(Cambiar esto para Carbonara)Chimera binaries could be installed automatically with the plugin after accepting ChimeraX licence terms,
but you can also link an existing installation. Default installation path assumed is *software/em/chimerax-1.0,
if you want to change it, set *CHIMERA_HOME* in *scipion.conf* file to the folder where ChimeraX is installed
or link your chimerax folder to *software/em/chimerax-1.0*.

- **Tests**

Tested with Carbonara version: XXXXX

To check the installation, simply run the following Scipion tests: 

- **Supported versions of CARBonAra**

XX,


=========
Protocols
=========

========
Examples
========

===
FAQ
===

===============
Buildbot status
===============

Status devel version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/carbonara_devel.svg

..
    Status production version: 

.. 
    image:: http://scipion-test.cnb.csic.es:9980/badges/carbonara_prod.svg





