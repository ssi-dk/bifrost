
## Workflow

Initialization -> Run

Initialization is the process of generating the appropriate run, sample, component, run_component, and sample_component structure for a "Run" of data. No data is actually processed at this time but the database is set up accordingly for the inputs.

A run in this system can be considered the same as a sequencing run or project. It's an organization structure for samples including which samples belong to the run and what components are planned to be run. Typically sequencing runs will form the basis of your data structure for the platform.

A sample in this system can be thought of as a single isolate from a single run. This platform is very sample orientated. Information related to the sample include read information, what components the sample will be run against, some properties such as species, and optional metadata. The idea is from a sample alone you can track all the components run against it and which ones it hasn't.

A component in this system can be thought of as a software pipeline. Each version of a component is considered a different component. There are several types of components, 

- Pipeline: The most standard which indicates running a software pipeline on a sample and record the results. Each pipeline should have a datadump file which summarizes output for the database.
- Stamper: A component which works exclusively off of the database. The idea behind this is if we need to create or adjust QC standards we can apply it quickly to all samples in the system without using the linux server to run new software.
- Connector: A connector converts the output of one to multiple components to fit a standard form. This can be useful for running visualization which requires set input on multiple pipelines. This one is currently not implemented.

On top of this components currently target 2 broad categories:

- Run (meant to affect whole runs and not individual samples)
- Sample (meant to affect individual samples)

A sample_component in this system is where you store the results of the specific sample and the specific component.

A run_component in this system is where you store the results of a specific run and a specific component.

Currently the bifrost system is only designed to handle Illumina Paired-end reads for bacterial WGS most nominally for QC related tasks.

## Starting up a run

01. find all potential samples which will be included in the run
02. create the samples folder based on the potential samples
03. find all potential components which will be included in the run
04. create the components folder based on the potential components
05. create the components into the component DB entry
06. create the run folder
07. save the samples and components (both sample and run) into the run DB entry
08. add run_components to run folder
09. add run_components DB entry with run and component info
10. add sample DB entry
11. add components to samples DB entry
12. add sample_components to sample folder
13. add sample_components DB entry with sample and component info
14. Create shell script for running program
15. run program

### Species picker logic

Possibilities:
  0. no provided species, no detected species -> leave blank
  1. no provided species -> set species to detected
  2. provided species not in DB -> set species to provided species
  3. provided species in DB as species -> set species to provided species
  4. provided species in DB as group
    a. detected species is a member of the group -> set species to detected species
    b. detected species is NOT a member of the group -> set species to None

## Anonymizer approach

New DB table with mapping from ID to name, names are internally saved as ID both in DB and in server a mirrored directory could exist linking the IDs to dynamix names which the DB mapping table would provide, updating the table should update the name in all locations then. If you don't have a mapping a name would be assigned. If you want to export data you could choose which to anonymize and which to leave.