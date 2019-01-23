import datetime
MINSPECIES = 0.95

MAXUNCLASSIFIED = 0.2

#CONTIGMAX = 1500

AVGCOVERAGE_FAIL = 10
AVGCOVERAGE_LOW = 25
AVGCOVERAGE_WARN = 50

NUMREADS_FAIL = 10000
#NUMREADS_WARN = 100000

MAXDIFFSIZE1x10x = 250000


def generate_summary(results, actions):
    summary = {}
    for result in results:
        summary[result["name"]] = result["status"] \
            + ":" + result["reason"] + ":" + str(result["value"])
    for component, action in actions.items():
        summary[component + ":action"] = action
    return summary


def evaluate_tests(tests, component):

    ## Check qcquickie results and define action
    results_dict = {result["name"]: result for result in tests}
    core_facility = False  # important
    supplying_lab = False
    for result in tests:
        if result["name"].startswith(component) or result["name"].startswith("whats_my_species"):
            if result["status"] == "fail" or result["status"] == "undefined":
                if result["effect"] == "supplying lab":
                    supplying_lab = True
                elif result["effect"] == "core facility":
                    core_facility = True
    if (component == "qcquickie" or component == "assemblatron") and \
    "whats_my_species:submitted == detected" in results_dict:
        if (results_dict[component + ":avgcoverage"]["status"] == "fail" and
            results_dict[component + ":avgcoverage"]["effect"] == "supplying lab") and \
            results_dict[component + ":1x10xsizediff"]["status"] == "fail" and \
            (results_dict["whats_my_species:detectedspeciesmismatch"]["status"] == "pass" and
             results_dict[component + ":1xgenomesize"]["status"] == "pass"):
                core_facility = True

    action = "OK"
    if supplying_lab:
        action = "supplying lab"
    if core_facility:
        action = "core facility"
    return action

def generate_stamp(actions, assemblatron):
    action = actions.get("assemblatron", "not checked")
    if action == "OK":
        value = "pass:OK"
    elif action == "not checked":
        value = "not checked:not checked"
    else:
        value = "fail:" + action
        
    return {
        "name": "ssi_stamper",
        "value": value,
        "date": datetime.datetime.utcnow()
    }

def test(whats_my_species, qcquickie, assemblatron, species, sample):


    results = []
    actions = {}

    try:
        test = {
            "name": "base:readspresent",
            "display_name": "No reads",
            "effect": "core facility",
            "value": ""
        }
        read_path = sample["reads"]["R1"]
        test["value"] = read_path
        if read_path == "":
            test["status"] = "fail"
            test["reason"] = "Read path is empty"
        else:
            test["status"] = "pass"
            test["reason"] = ""
    except KeyError:
        test["status"] = "fail"  # instead of KeyError
        test["reason"] = "Read path is undefined (KeyError)"
    results.append(test)

    if whats_my_species:
        try:
            test = {
                "name": "whats_my_species:minspecies",
                "display_name": "Multiple species detected",
                "effect": "supplying lab",
                "value": ""
            }
            detected = whats_my_species["summary"]["percent_classified_species_1"] + \
                whats_my_species["summary"]["percent_unclassified"]
            test["value"] = detected  # Add percentage
            if isinstance(test["value"], float):
                test["value"] = round(test["value"], 3)
            if detected < MINSPECIES:
                test["status"] = "fail"
                test["reason"] = "Value ({}) is below threshold ({})".format(
                    test["value"], MINSPECIES)
            else:
                test["status"] = "pass"
                test["reason"] = ""
        except KeyError:
            test["status"] = "fail"
            test["reason"] = "Detected value is undefined (KeyError)"
        results.append(test)

        try:
            test = {
                "name": "whats_my_species:maxunclassified",
                "display_name": "High unclassified",
                "effect": "supplying lab",
                "value": ""
            }
            unclassified = whats_my_species["summary"]["percent_unclassified"]
            test["value"] = unclassified
            if isinstance(test["value"], float):
                test["value"] = round(test["value"], 3)
            if unclassified >= MAXUNCLASSIFIED:
                test["status"] = "fail"
                test["reason"] = "Value ({}) is above threshold ({})".format(
                    test["value"], MAXUNCLASSIFIED)
            else:
                test["status"] = "pass"
                test["reason"] = ""
        except KeyError:
            test["status"] = "fail"
            test["reason"] = "Unclassified value is undefined (KeyError)"
        results.append(test)


        # This one is not in the reporter as it shouldnt appear often and there is no
        # obvious place to put it
        try:
            test = {
                "name": "whats_my_species:nosubmitted",
                "display_name": "No species submitted - using default values",
                "effect": "supplying lab",
                "value": ""
            }
            detected = species["organism"]
            test["value"] = detected
            if detected == "default":
                test["status"] = "fail"
                test["reason"] = "Detected species not in bifrost db. Please report this to system admin."
            else:
                test["status"] = "pass"
                test["reason"] = ""
        except KeyError:
            test["status"] = "fail"
            test["reason"] = "Detected species is undefined (KeyError)"
        results.append(test)

        # Submitted == detected
        try:
            test = {
                "name": "whats_my_species:detectedspeciesmismatch",
                "display_name": "Detected species mismatch",
                "effect": "supplying lab",
                "value": ""
            }
            detected = whats_my_species["summary"]["name_classified_species_1"]
            submitted = sample["properties"]["provided_species"]
            test["value"] = "{}=={}".format(submitted, detected)

            if submitted is None:
                test["status"] = "fail"
                test["reason"] = "No submitted species"
            else:
                if submitted != detected and submitted != species["group"]: # If submitted is group it should be contained 
                    test["status"] = "fail"
                    test["reason"] = "Detected species ({}) different than expected ({})".format(detected, submitted)
                else:
                    test["status"] = "pass"
                    test["reason"] = ""
        except KeyError as e:
            test["status"] = "fail"
            test["reason"] = "Value is undefined (KeyError): {}".format(
                e.args[0])
        results.append(test)

    for assembly_component, comp_name in ((qcquickie, "qcquickie"), (assemblatron, "assemblatron")):
        if assembly_component and len(assembly_component):
            # Genome size check for 1x
            try:
                test = {
                    "name": comp_name + ":1xgenomesize",
                    "display_name": "Atypical genome size (1x)",
                    "effect": "supplying lab",
                    "value": ""
                }
                size = int(assembly_component["summary"]["bin_length_at_1x"])
                min_size = int(species["min_length"])
                max_size = int(species["max_length"])
                test["value"] = size

                if min_size < size < max_size:
                    test["status"] = "pass"
                    test["reason"] = ""
                else:
                    test["status"] = "fail"
                    test["reason"] = "Value ({}) below or above expected ({}, {})".format(
                        size, min_size, max_size)
            except KeyError:
                test["status"] = "fail"
                test["reason"] = "1x Genome Size is undefined (KeyError)"
            results.append(test)

            # Genome size check for 10x
            try:
                test = {
                    "name": comp_name + ":10xgenomesize",
                    "display_name": "Atypical genome size (10x)",
                    "effect": "supplying lab",
                    "value": ""
                }
                size = int(assembly_component["summary"]["bin_length_at_10x"])
                min_size = int(species["min_length"])
                max_size = int(species["max_length"])
                test["value"] = size

                if min_size < size < max_size:
                    test["status"] = "pass"
                    test["reason"] = ""
                else:
                    test["status"] = "fail"
                    test["reason"] = "Value ({}) below or above expected ({}, {})".format(
                        size, min_size, max_size)
            except KeyError:
                test["status"] = "fail"
                test["reason"] = "10x Genome Size is undefined (KeyError)"
            results.append(test)

            # Genome size difference 1x 10x
            try:
                test = {
                    "name": comp_name + ":1x10xsizediff",
                    "display_name": "Atypical genome size difference (1x - 10x)",
                    "effect": "supplying lab",
                    "value": ""
                }
                size_1x = assembly_component["summary"]["bin_length_at_1x"]
                size_10x = assembly_component["summary"]["bin_length_at_10x"]
                diff = size_1x - size_10x
                max_diff = MAXDIFFSIZE1x10x
                if isinstance(diff, float):
                    test["value"] = "{}".format(round(diff, 3))
                else:
                    test["value"] = "{}".format(diff)
                if diff < max_diff:
                    test["status"] = "pass"
                    test["reason"] = ""
                else:
                    test["status"] = "fail"
                    test["reason"] = "Value ({}) above expected ({})".format(
                        diff, max_diff)
            except KeyError:
                test["status"] = "fail"
                test["reason"] = "Value is undefined (KeyError)"
            results.append(test)

            # Average coverage
            try:
                test = {
                    "name": comp_name + ":avgcoverage",
                    "display_name": "Coverage test",
                    "effect": "core facility",
                    "value": ""
                }
                coverage = assembly_component["summary"]["bin_coverage_at_1x"]
                test["value"] = coverage
                test["value"] = round(test["value"], 3)

                if coverage < AVGCOVERAGE_FAIL:
                    test["status"] = "fail"
                    test["reason"] = "Lack of reads ({} < {})".format(
                        coverage, AVGCOVERAGE_FAIL)
                elif coverage < AVGCOVERAGE_LOW:
                    test["status"] = "fail"
                    test["reason"] = "Not enough reads ({} < {})".format(
                        coverage, AVGCOVERAGE_LOW)
                elif coverage < AVGCOVERAGE_WARN:
                    test["status"] = "fail"
                    test["reason"] = "Low reads ({} < {})".format(
                        coverage, AVGCOVERAGE_WARN)
                    test["effect"] = "supplying lab"
                else:
                    test["status"] = "pass"
                    test["reason"] = ""
            except KeyError:
                test["status"] = "fail"
                test["reason"] = "Coverage is undefined (KeyError)"
            results.append(test)

            # Minimum number of reads
            try:
                test = {
                    "name": comp_name + ":numreads",
                    "display_name": "Number of filtered reads below minimum",
                    "effect": "core facility",
                    "value": ""
                }
                num_reads = assembly_component["summary"]["filtered_reads_num"]
                test["value"] = assembly_component["summary"]["filtered_reads_num"]
                if num_reads < NUMREADS_FAIL:
                    test["status"] = "fail"
                    test["reason"] = "Filtered reads below minimum ({} < {})".format(
                        num_reads, NUMREADS_FAIL)
                else:
                    test["status"] = "pass"
                    test["reason"] = ""
            except KeyError:
                test["status"] = "fail"
                test["reason"] = "Number of reads is undefined (KeyError)"
            results.append(test)
            actions[comp_name] = evaluate_tests(
                results, comp_name)

    return results, generate_summary(results, actions), generate_stamp(actions, assemblatron)


