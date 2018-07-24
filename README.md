SerumQC-private

## Running SerumQC:

Activate conda environment:
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

Set the components you want to run (component names are the filenames in /snakefiles dir), default is `qcquickie,assembly,analysis`.

Use `use_mongodb=False` in config if you don't want to store the run, samples and results in the database.

```
snakemake -s code/Snakefile --config components=qcquickie
```
