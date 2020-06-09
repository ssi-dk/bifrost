#!/usr/bin/env python3
"""
Launcher file for accessing dockerfile commands
"""
import argparse
import json
import subprocess
import os
import sys
import traceback
from bifrostlib import datahandling

COMPONENT: dict = datahandling.load_yaml(os.path.join(os.path.dirname(__file__), 'config.yaml'))


def parse_args():
    """
    Arg parsing via argparse
    """
    description: str = (
        f"-Description------------------------------------\n"
        f"{COMPONENT['details']['description']}"
        f"------------------------------------------------\n\n"
        f"*Run command************************************\n"
        f"docker run \ \n"
        f" -e BIFROST_DB_KEY=mongodb://<user>:<password>@<server>:<port>/<db_name> \ \n"
        f" -v <input_path>:/input \ \n"
        f" -v <output_path>:/output \ \n"
        f" {COMPONENT['dockerfile']} \ \n"
        f"************************************************\n"
    )
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--install',
                        action='store_true',
                        help='Install component')
    parser.add_argument('-info', '--info',
                        action='store_true',
                        help='Provides basic information on component')
    parser.add_argument('-pre', '--pre_script',
                        help='Pre script template run before sample script')
    parser.add_argument('-per', '--per_sample_script',
                        help='Per sample script template run on each sample')
    parser.add_argument('-post', '--post_sample_script',
                        help='Post script template run after sample script')
    parser.add_argument('-meta', '--run_metadata',
                        required=True,
                        help='Run metadata tsv')
    parser.add_argument('-metamap', '--run_metadata_column_remap',
                        help='Remaps metadata tsv columns to bifrost values')

    args: argparse.Namespace = parser.parse_args()

    if not datahandling.check_db_connection_exists():
        message: str = (
            f"ERROR: Connection to DB not establised.\n"
            f"please ensure env variable BIFROST_DB_KEY is set and set properly\n"
        )
        print(message)
    else:
        print(datahandling.get_connection_info())

    if args.info:
        show_info()
    elif args.install:
        install_component()
    elif args.sample_id is not None:
        run_pipeline(args)


def show_info():
    """
    Shows information about the component
    """
    message: str = (
        f"Component: {COMPONENT['name']}\n"
        f"Version: {COMPONENT['version']}\n"
        f"Details: {json.dumps(COMPONENT['details'], indent=4)}\n"
        f"Requirements: {json.dumps(COMPONENT['requirements'], indent=4)}\n"
        f"Output files: {json.dumps(COMPONENT['db_values_changes']['files'], indent=4)}\n"
    )
    print(message)


def install_component():
    component: list[dict] = datahandling.get_components(component_names=[COMPONENT['name']], component_versions=[COMPONENT['version']])
    if len(component) == 1:
        print(f"Component has already been installed")
    elif len(component) > 1:
        print(f"Component exists multiple times in DB, please contact an admin to fix this in order to proceed")
    else:
        datahandling.post_component(COMPONENT)
        component: list[dict] = datahandling.get_components(component_names=[COMPONENT['name']], component_versions=[COMPONENT['version']])
        if len(component) != 1:
            print(f"Error with installation of {COMPONENT['name']} v:{COMPONENT['version']} \n")
            exit()


def run_pipeline(args: object):
    """
    Runs pipeline
    """
    sample: list[dict] = datahandling.get_samples(sample_ids=[args.sample_id])
    component: list[dict] = datahandling.get_components(component_names=[COMPONENT['name']], component_versions=[COMPONENT['version']])
    if len(component) == 0:
        print(f"component not found in DB, installing it:")
        datahandling.post_component(COMPONENT)
        component: list[dict] = datahandling.get_components(component_names=[COMPONENT['name']], component_versions=[COMPONENT['version']])
        if len(component) != 1:
            print(f"Error with installation of {COMPONENT['name']} v:{COMPONENT['version']} \n")
            exit()

    else:
        print(f"snakemake -s /bifrost/{COMPONENT['name']}/pipeline.smk --config sample_id={str(sample[0]['_id'])} component_id={str(component[0]['_id'])}")
        try:
            process: subprocess.Popen = subprocess.Popen(
                f"pipeline.smk sample_id={str(sample[0]['_id'])} component_id={str(component[0]['_id'])}",
                stdout=sys.stdout,
                stderr=sys.stderr,
                shell=True
            )
            process.communicate()
        except:
            print(traceback.format_exc())


if __name__ == '__main__':
    parse_args()
