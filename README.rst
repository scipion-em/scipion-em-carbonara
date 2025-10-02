=======================
Carbonara scipion plugin
=======================

This plugin allows to use run **CARBonAra** commands within the **Scipion** framework.

CARBonAra allows protein sequence designing based on structure atomic coordinates. Since this method takes into account the context in sequence generation, it is particularly helpful to design proteins embedded in molecular environments that include non-protein entities such as nucleic acids (See `<https://github.com/LBM-EPFL/CARBonAra/>`_  for details and features).

==========================
Download the repository
==========================
.. code-block::

    git clone git@github.com:scipion-em/scipion-em-carbonara.git

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
      
OR through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version**
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





