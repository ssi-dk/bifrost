import stamps.ssi_stamp
import traceback
from bifrostlib import datahandling


def get_components(sample_dict):
    component_names = ["whats_my_species", "qcquickie", "assemblatron"]
    comps = {}
    for component in component_names:
        comps[component] = datahandling.get_sample_components(
            sample_ids=[str(sample_dict["_id"])],
            component_names=[component])[0]
    return comps


def script__test_ssi_stamper(sample, sample_yaml, sample_component, log_err):
    # Genering error handling to redirect output to stderr file
    try:
        comps = get_components(sample)
        if "detected_species" in sample["properties"]:
            species = datahandling.load_species(sample["properties"]["detected_species"])
        else:
            species = datahandling.load_species(None)
        
        results, summary, stamp = stamps.ssi_stamp.test(comps['whats_my_species'], comps['qcquickie'], comps['assemblatron'],
                                species, sample)

        datadump_dict = datahandling.load_sample_component(sample_component)
        datadump_dict["summary"] = summary
        datadump_dict["results"] = results
        datahandling.save_sample_component(datadump_dict, sample_component)

        # Get the _id back

        datadump_dict = datahandling.load_sample_component(sample_component)
        stamp["_sample_component"] = datadump_dict["_id"]

        stamp_dict = sample.get("stamps", {})
        stamp_list = stamp_dict.get("stamp_list", [])
        stamp_list.append(stamp)
        stamp_dict["stamp_list"] = stamp_list
        stamp_dict[stamp["name"]] = stamp
        sample["stamps"] = stamp_dict

        datahandling.save_sample(sample, sample_yaml)
        return 0
    except Exception:
        datahandling.log(log_err, str(traceback.format_exc()))
        exit(1)

script__test_ssi_stamper(snakemake.params.sample,
                        snakemake.params.sample_yaml,
                        snakemake.params.sample_component,
                        snakemake.log.err_file)
