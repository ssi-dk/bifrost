import re
import pandas
from ruamel.yaml import YAML
import sys

# should scan samples.yaml and determine what I can report on

configfile: os.path.join(os.path.dirname(workflow.snakefile), "../config.yaml")
# requires --config R1_reads={read_location},R2_reads={read_location}
sample = config["Sample"]
global_threads = config["threads"]
global_memory_in_GB = config["memory"]

yaml = YAML(typ='safe')
yaml.default_flow_style = False

with open(sample, "r") as sample_yaml:
    config_sample = yaml.load(sample_yaml)


# run minireport