import sys
import argparse
from bifrostlib import datahandling

# Argument parsing
parser = argparse.ArgumentParser(
    description=("This script checks the status for all the sample components "
                 "in a run. It is used against a test run to make sure the"
                 " tested version of bifrost works. Any status other than"
                 " 'Success' or 'requirements not met' in the expected counts "
                 "will fail the test. This uses bifrostlib and checks the "
                 "database following the key in the env var BIFROST_DB_KEY")
)
parser.add_argument("run_name")
parser.add_argument("expected_success_count",
    help="Expected number of sample_components with 'Success' status",
    type=int)
parser.add_argument("expected_req_count",
    help="Expected number of sample_components with 'requirements not met' status",
    type=int)

args = parser.parse_args()

# Testing assumptions

run = datahandling.get_runs(names=[args.run_name])
if not len(run):
    sys.stderr.write("Test failed. Run {} not found.\n".format(args.run_name))
    exit()

run = run[0]

sample_ids = [str(s["_id"]) for s in run["samples"]]

s_cs = datahandling.get_sample_components(sample_ids=sample_ids)

if not len(s_cs):
    sys.stderr.write(
        "Test failed. No sample_components found.\n")
    exit()

status_count = {}
for s_c in s_cs:
    status = s_c["status"]
    status_count[status] = status_count.get(status, 0) + 1

fail = False
for key, val in status_count.items():
    if key == "Success" and val == args.expected_success_count:
        continue
    elif key == "requirements not met" and val == args.expected_req_count:
        continue
    else:
        fail = True
        break

if fail:
    sys.stderr.write(
        ("Test failed. Expectations:\nSuccess: {}\nrequirements not met:"
         " {}\nFound:\n{}\n").format(
             args.expected_success_count,
             args.expected_req_count,
             status_count
         ))
    exit()
else:
    sys.stdout.write("Test passed.\n")



