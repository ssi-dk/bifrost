# Welcome to bifrost's documentation

Bifrost is a bacterial quality control (QC) software platform developed by Statens Serum Institut (SSI) for tracking collections of samples at the sample level for detail. By tracking what analyses have been run against each sample in a centralized database (DB) we can perform rapid routine analysis against all values that enter the system. To ensure that the analyses are accurate it is important to run and track each sample through QC which is appropriate for the downstream analysis you want to perform and is thus a central component of bifrost. A centralized DB also gives greater power in analyzing your sample collection through either direct DB interaction or via bifrost's web based graphical user interface (gui). Bifrost provides a platform allowing for real-time information on how a new sample compares against all existing DB samples enabling more data-based approaches for analysis. Bifrost has been implemented at SSI and we are as of this writing in the process to rolling this platform out to all microbiology labs in Denmark and will see continued development to fulfill the needs of software as infrastructure.

## High Level Concepts

Bifrost is built around the idea of tracking analysis of samples through a DB with a sample being the basepoint. This means

## DB Structure

The structure of bifrost is currently based around the following tables:

- run
- - Here we have a continuous flow of samples that are to be used for surveillance which undergo WGS. This is the basis for a run structure. A run structure can also be used to describe projects
- sample
- component
- run_component
- sample_component

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

## Installing bifrost

### Requirements

Python >=3.7

anaconda/miniconda >=4.3

MongoDB >=4.0.0

pip >= 9.0

### Setting up MongoDB

Get MongoDB (community edition) and set it up following their documentation. 
You'll need a database and read/write account for it.

https://www.mongodb.com/download-center/community

The system relies on a species collection containing species-specific information used for 
QC testing and species identification. 

Here is an example of an entry:

```
{
    "_id" : ___,
    "organism" : "E. coli",
    "group" : "Escherichia",
    "mlst_species" : "/path_to_db/ariba/mlst/Escherichia_coli_1/ref_db/",
    "ncbi_species" : "Escherichia coli",
    "min_length" : 4500000, #min and max genome size to pass QC
    "max_length" : 5800000
}
```

After this step you should have a mongodb instance running, a user and a species collection with the species you expect to find in the samples.

### Setting up the conda environment

Install the conda environment using envs/bifrost_for_install.yaml

### Installing bifrostlib

Install bifrostlib withing the conda environment:

```
source activate bifrost
pip install lib/bifrostlib/
```


### Adjusting config

Adjust the main config.yaml and the individual components config.yaml.

### Setting up the database key file location environment variable

Create a file named `keys.txt` with the mongodb URI including the database to use.

Create an environment variable with the mongodb key location.
```
export BIFROST_DB_KEY=/path/to/keys.txt
```

You can also set the key on the conda environment so it gets set when the environment is active.

https://conda.io/docs/user-guide/tasks/manage-environments.html#saving-environment-variables

## Running bifrost:

Activate conda environment `bifrost`:
```
source activate bifrost
```

Create output directory, and move into it 
```
mkdir run_name # Run name will be taken from the directory name or it can be set through config
cd run_name
```

Link the demultiplexed (fastq.gz) samples folder as samples (by default):
```
ln -s path/to/run samples
```

Copy/clone this directory in src directory:

```
git clone git@github.com:ssi-dk/bifrost.git
mv bifrost src
```

Run snakemake to generate the run command and the folder structures.

Set the components you want to run (component names are the filenames in /components dir), default is `qcquickie,assemblatron,analyzer,ssi_stamper`.

Use `use_mongodb=False` in config if you don't want to store the run, samples and results in the database.
For now however, you need to use the database to access the species table to run mlst.

Then start the program running:

```
snakemake -s src/bifrost.smk
```

You can change any option in config.yaml using the `--config` argument. Here is an example with some useful ones:

```
--config components=qcquickie,assemblatron,analyzer partition=project run_name=myrun assembly_with=spades
```

For more info on each parameter, check the sample config file.
