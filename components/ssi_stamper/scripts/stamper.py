import sys
import os
import datetime
import argparse
from bifrostlib import datahandling

from stamps import ssi_stamp


def run_test(sample_dict, stamp):
    if stamp == "ssi_stamper":
        component_names = ["whats_my_species", "qcquickie", "assemblatron"]
        comps = {}
        for component in component_names:
            comps[component] = datahandling.load_last_sample_component(
                str(sample_dict["_id"]), component)
        species = datahandling.load_species(sample_dict["properties"].get("species", None))
        return ssi_stamp.test(comps["whats_my_species"], comps["qcquickie"], comps["assemblatron"],
                            species, sample_dict)


def get_git_hash():
    head_path = os.path.join(os.path.dirname(__file__), "../.git/HEAD")
    with open(head_path) as head:
        head_line = head.readline().strip()
        loc = head_line.split(" ")[1]
    hash_path = os.path.join(os.path.dirname(__file__), "../.git/" + loc)
    with open(hash_path) as g_hash_file:
        return g_hash_file.readline().strip()

def create_component(stamp_name):
    return {
        "name": stamp_name,
        "version": "1.0",
        "git": get_git_hash()
    }

def create_sample_component(sample, component, results, summary):
    return {
        "sample": {
            "_id": sample["_id"],
            "name": sample["name"]
        },
        "component": {
            "_id": component["_id"],
            "name": component["name"]
        },
        "results": results,
        "summary": summary,
        "status": "Success",
        "setup_date": datetime.datetime.utcnow()
    }

def main(argv):

    parser = argparse.ArgumentParser(description='Run a stamp against samples in the DB.')
    parser.add_argument('stamp', choices=["ssi_stamper"],
                        help='name of the stamp to use')
    parser.add_argument('-s', '--samples', 
                        help='list of comma-separated sample ids. Use this OR --runs OR --all')
    parser.add_argument('-r', '--runs',
                        help='list of comma-separated run ids. Use this OR --samples OR --all')
    parser.add_argument('-A', '--all',
                        help='Test all samples. Use this OR --samples OR --runs', action = "store_true")

    args = parser.parse_args()

    stamp_name = args.stamp
    more_than_one_arg = 0

    if args.all:
        more_than_one_arg += 1
        print("This will apply the stamp on all the samples.")
        if input("Are you sure? ") == "y":
            sample_ids = list(map(str, datahandling.load_all_samples()))
        else:
            exit()
    if args.samples:
        more_than_one_arg += 1
        sample_ids = args.samples.split(",")
    if args.runs:
        more_than_one_arg += 1
        run_ids = args.runs.split(',')
        sample_ids = datahandling.load_samples_from_runs(run_ids)
        sample_ids = list(map(str, sample_ids))
    
    if more_than_one_arg > 1:
        exit("Invalid parameters. Use --samples OR --runs OR all")

    print("Testing samples ({}):\n{}".format(len(sample_ids), "\n".join(sample_ids)))

    if len(sample_ids) == 0:
        exit()
    
    component_db = datahandling.save_component_to_db(
        create_component(stamp_name))

    for sample_id in sample_ids:
        print("Testing sample: " + sample_id)
        sample_db = datahandling.load_sample_from_db(sample_id)
        results, summary, stamp = run_test(sample_db, stamp_name)

        print("Tested sample {}.\n{}\n\n".format(sample_id, stamp))
        # Saving stuff
        
        s_c = create_sample_component(sample_db, component_db, results, summary)
        s_c_db = datahandling.save_sample_component_to_db(s_c)

        if "components" in sample_db:
            sample_db["components"].append({"name": component_db["name"], "_id": component_db["_id"]})
        else:
            sample_db["components"] = [{"name":component_db["name"], "_id":component_db["_id"]}]
        
        stamp["_sample_component"] = s_c_db["_id"]
        if "stamps" in sample_db:
            sample_db["stamps"].append(stamp)
        else:
            sample_db["stamps"] = [stamp]

        datahandling.save_sample_to_db(sample_db)

if __name__ == "__main__":
    main(sys.argv)
