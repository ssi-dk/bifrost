.. bifrost documentation master file, created by
   sphinx-quickstart on Thu Dec 13 14:56:13 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================================
Welcome to bifrost's documentation!
===================================

Bifrost is a bacterial quality control (QC) software platform developed by Statens Serum Institut (SSI) for tracking collections of samples at the sample level for detail. By tracking what analyses have been run against each sample in a centralized database (DB) we can perform rapid routine analysis against all values that enter the system. To ensure that the analyses are accurate it is important to run and track each sample through QC which is appropriate for the downstream analysis you want to perform and is thus a central component of bifrost. A centralized DB also gives greater power in analyzing your sample collection through either direct DB interaction or via bifrost's web based graphical user interface (gui). Bifrost provides a platform allowing for real-time information on how a new sample compares against all existing DB samples enabling more data-based approaches for analysis. Bifrost has been implemented at SSI and we are as of this writing in the process to rolling this platform out to all microbiology labs in Denmark and will see continued development to fulfill the needs of software as infrastructure.

===================
High level concepts
===================

Bifrost is built around the idea of tracking analysis of samples through a DB with a sample being the basepoint. This means

============
DB Structure
============

The structure of bifrost is currently based around the following tables
* run
    * Here we have a continuous flow of samples that are to be used for surveillance which undergo WGS. This is the basis for a run structure. A run structure can also be used to describe projects
* sample
* component
* run_component
* sample_component

=================
Starting up a run
=================


# species picker logic
# Possibilities, 
#   0. no provided species, no detected species -> leave blank
#   1. no provided species -> set species to detected
#   2. provided species not in DB -> set species to provided species
#   3. provided species in DB as species -> set species to provided species
#   4. provided species in DB as group
#       a. detected species is a member of the group -> set species to detected species
#       b. detected species is NOT a member of the group -> set species to None

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
