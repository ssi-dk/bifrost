SerumQC-private

## Running SerumQC:

Activate conda environment. `env_qcquickie` is the name of the environment in the SSI servers:
```
source activate env_qcquickie
```

Create output directory, and move into it 
```
mkdir run_name # Or whatever folder name you prefer
cd run_name
```

Link the demultiplexed samples folder as samples (by default):
```
ln -s path/to/run samples
```

Copy/clone this directory in code:

```
cp -r path/to/this/repo code
```

Run snakemake to generate the run command and the folder structures.

Set the components you want to run (component names are the filenames in /snakefiles dir), default is `qcquickie,assemblatron,analyzer`.

Use `use_mongodb=False` in config if you don't want to store the run, samples and results in the database.
For now however, you need to use the database to access the species table to run mlst.
Make sure you have the mongo url in the file: `/resources/keys.txt`

Then start the program running:

```
snakemake -s code/Snakefile --config components=qcquickie
```

You can change any option in config.yaml using the `--config` argument. Here is an example with some useful ones:

```
--config components=qcquickie,assemblatron,analyzer partition=project run_name=myrun assembly_with=spades
```

For more info on each parameter, check the sample config file.
