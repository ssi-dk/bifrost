# Bifrost KMA README

This is the basic folder structure with a sample run containing only one sample:

```
├── config
│   ├── config.yaml
│   └── keys.txt
├── input
│   └── 2019
│       └── sample_run
│           ├── SAMPLE1_S2_L001_R1_001.fastq.gz
│           └── SAMPLE1_S2_L001_R2_001.fastq.gz
├── output
│   └── 2019
│       └── sample_run
│           ├── bifrost
│           ├── components
│           ├── config.yaml
│           ├── SAMPLE1
│           ├── run_cmd_bifrost.sh
│           ├── samples
│           └── src
└── README.md
```

`config` contains `config.yaml` which includes KMA specific config values to run the pipeline such as the group name for the torque queue and which components to run.

## Basic workflow

The basic workflow starts when the user adds a new directory (considered a 'run') to the input/year directory.
This run directory should contain the read pairs with names matching the regex specified in the config file. The default name from illumina MySeq machines works.

The next step is either triggered by the user or, when implemented, executed by a chron job.

## Setting up the run

- Create a directory in output/year (we'll call it output/year/run, an example could be output/2019/RUN2019001) with the same name as the one created in input/year (input/year/run)
- Symlink input/year/run into output/year/run/samples
- Copy recursively the source code located in /home/projects/ssi_10003/apps/bifrost-private as output/year/run/src
- If you have a sample sheet with sample metadata, add it as output/year/run/sample_sheet.xlsx (or sample_sheet.tsv, but you need to set config to sample_sheet=sample_sheet.tsv)

### Sample sheet





## Running bifrost

- Load the environment:

```
module load tools
module load anaconda3/4.4.0
source activate bifrost
```

- Set the location of the MongoDB key file. The file should be named keys.txt and contain the mongo URI for the KMA

```
export BIFROST_DB_KEY=/home/projects/KMA_project/data/bifrost/config/keys.txt
```

- Run bifrost from /output/year/run with the following command:

```
snakemake -s src/bifrost.smk --configfile /home/projects/KMA_project/data/bifrost/config/config.yaml
```

You can check the running status in your run checker.

The output for each sample will be generated in output/year/run/samplename 

## Debugging

If there are any problems running the pipeline contact Martin (mbas@ssi.dk) or Kim (kimn@ssi.dk). However, you can access some logs if you want to get more detail on where the program is failing.

There are two main stages in the execution: Creating the run and running the components.

### Creating the run

If you can't find the run in the run checker, bifrost probably hasn't finished creating it. Another way to confirm is to go to a sample folder (if they exist) and check that there are no component folders created yet.

The logs for this stage can be found in output/year/run/.snamemake/log/ and in output/year/run/bifrost/log

### Running the components

If a component fails, (it shows as Failed or in rare cases as Running for a long time, in the Run checker), you can check the logs in output/year/run/<sample>/<component>/log/

There might be some logs from torque in your home directory, although we plan to change it to a proper location.
