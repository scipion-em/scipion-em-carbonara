=========================
CARBonAra Scipion plugin
=========================

This plugin allows to run **CARBonAra** commands within the **Scipion** framework.

CARBonAra allows protein sequence designing based on structure atomic coordinates.
Since this method takes into account the context to perform sequence generation, it is particularly helpful to design proteins
embedded in molecular environments that include non-protein entities such as nucleic acids
(See `<https://github.com/LBM-EPFL/CARBonAra/>`_  for details and features).

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

      cd scipion-em-carbonara
      scipion3 installp -p scipion-em-carbonara
      
OR through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**
      
- **Developer's version**

.. code-block:: 

      cd scipion-em-carbonara
      scipion3 installp -p path_to_scipion-em-carbonara --devel

- **Configuration variables**

The installation of CARBonAra within Scipion automatically defines the variable `CARBONARA_ENV_ACTIVATION` for specifying
how to activate the conda environment required to run CARBonAra and wil be added to the final carbonara command. 

.. code-block::

    CARBONARA_ENV_ACTIVATION = conda activate carbonara
    
- **Binary files**

The carbonara executable for Scipion will be automatically installed in your folder miniconda/envs/carbonara/bin.

======
Tests
======
Tested with the first published Carbonara version. To check the installation, simply run the following Scipion test:

.. code-block:: 

      scipion tests carbonara.tests.test_protocol_carbonara_sampling_sequence

==========
Protocols
==========
* protocol_sampling_sequence: Automatic method to generate multiple sequences able to fold in a certain structural configuration.

=========
Examples
=========

====
FAQ
====

================
Buildbot status
================

Status devel version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/carbonara_devel.svg

..
    Status production version: 

.. 
    image:: http://scipion-test.cnb.csic.es:9980/badges/carbonara_prod.svg





