# Bifrost KMA README

This is the basic folder structure with a sample run containing only one sample:

```bash
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
├── sample_sheet.tsv
└── README.md <-- this file
```

`config` contains `config.yaml` which includes KMA specific config values to run the pipeline such as the group name for the torque queue and which components to run.
The sample_sheet.tsv is a template you can use when setting up a run (see below) to add metadata to the sample.

## Before starting

Create a `.condarc` in your home directory with the following content so that the bifrost environment can be run:

```
envs_dirs:
  - /home/projects/ssi_10003/apps/anaconda/envs

pkgs_dirs:
  - /home/projects/ssi_10003/apps/anaconda/pkgs
```

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

You can use /sample_sheet.tsv as an example. It should have a header and the order is not important as it is being read by pandas. The fields are:

**SampleID:** Mandatory. This name is mapped in the database to sample_name, but our code currently depends on that name instead of the mapped one as an oversight. I'll add that to the task list so in the future it will be sample_name (or whatever field you configure in config.yaml).

**group:** Recommended. We use it for the different supplying labs and it can be used to sort samples in the web reporter, but it is not required if it doesn't really apply to your case.

**provided_species:** Strongly recommended. We detect the species anyway, but it is useful for QC. The current behavior is that qc will flag it if missing, or if it is different than the detected one (to identify contaminations or wrongly identified species). The QC behavior can be changed of course.

**Comments:** Optional. This value is shown in the reporter, so don't store anything here that you don't want anyone with access to the web reporter to see. In our case, it is empty in most samples but we use it for library preparation failures and misc information. As with SampleID, we will change the name to comments to have a consistent naming scheme.

**emails:** Optional. Shown as well in the reporter, we use it to show a contact person for the sample. An email string or a comma separated list of emails is the expected value.

You can add any extra values to the sample sheet and they will not be shown in the reporter, but will be stored in the database (in the sample object, under sample_sheet).

Also note, to make it work with the tsv you need to change the default config.yaml value for sample_sheet to `sample_sheet=sample_sheet.tsv`.

## Running bifrost

- Load the environment:

```bash
module load tools
module load anaconda3/4.4.0
source activate bifrost
```

- Set the location of the MongoDB key file. The file should be named keys.txt and contain the mongo URI for the KMA

```bash
export BIFROST_DB_KEY=/home/projects/KMA_project/data/bifrost/config/keys.txt
```

- Run bifrost from /output/year/run with the following command:

```bash
snakemake -s src/bifrost.smk --configfile /home/projects/KMA_project/data/bifrost/config/config.yaml
```

You can check the running status in your run checker.

The output for each sample will be generated in output/year/run/samplename

## Maintenance

### Removing a run

There is a script in scripts/db_management/remove_run.py. Make sure to have the bifrost environment enabled and your correct key path in the BIFROST_DB_KEY variable. Then run the script:

```bash
python remove_run.py 5b894291104fd46ecabe6638
```

Where 5b894291104fd46ecabe6638 is the id for your run (without the "ObjectId()" part). You can find the id of the run in the run folder/bifrost/run.yaml or in the database directly. To connect to the database, you need the mongo client which you can find in computerome already by running:

module load tools
module load mongodb/3.4.4

and then running mongo. You can use the credentials in your key file to log in and take a look.

## Debugging

If there are any problems running the pipeline contact Martin (mbas@ssi.dk) or Kim (kimn@ssi.dk). However, you can access some logs if you want to get more detail on where the program is failing.

There are two main stages in the execution: Creating the run and running the components.

### Creating the run

If you can't find the run in the run checker, bifrost probably hasn't finished creating it. Another way to confirm is to go to a sample folder (if they exist) and check that there are no component folders created yet.

The logs for this stage can be found in output/year/run/.snamemake/log/ and in output/year/run/bifrost/log

### Running the components

If a component fails, (it shows as Failed or in rare cases as Running for a long time, in the Run checker), you can check the logs in `output/year/run/<sample>/<component>/log/`

There might be some logs from torque in your home directory, although we plan to change it to a proper location.
