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

========
Workflow
========

Initialization -> Run
    Initialization is the process of generating the appropriate run, sample, component, run_component, and sample_component structure for a "Run" of data. No data is actually processed at this time but the database is set up accordingly for the inputs.

    A run in this system can be considered the same as a sequencing run or project. It's an organization structure for samples including which samples belong to the run and what components are planned to be run. Typically sequencing runs will form the basis of your data structure for the platform.

    A sample in this system can be thought of as a single isolate from a single run. This platform is very sample orientated. Information related to the sample include read information, what components the sample will be run against, some properties such as species, and optional metadata. The idea is from a sample alone you can track all the components run against it and which ones it hasn't.

    A component in this system can be thought of as a software pipeline. Each version of a component is considered a different component. There are several types of components, 
        Pipeline: The most standard which indicates running a software pipeline on a sample and record the results. Each pipeline should have a datadump file which summarizes output for the database. 
        Stamper: A component which works exclusively off of the database. The idea behind this is if we need to create or adjust QC standards we can apply it quickly to all samples in the system without using the linux server to run new software.
        Connector: A connector converts the output of one to multiple components to fit a standard form. This can be useful for running visualization which requires set input on multiple pipelines. This one is currently not implemented.

        On top of this components currently target 2 broad categories:
        Run (meant to affect whole runs and not individual samples)
        Sample (meant to affect individual samples)

    A sample_component in this system is where you store the results of the specific sample and the specific component.

    A run_component in this system is where you store the results of a specific run and a specific component.

Currently the bifrost system is only designed to handle Illumina Paired-end reads for bacterial WGS most nominally for QC related tasks.

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
