import stamps.ssi_stamp
from bifrostlib import datahandling


def get_components(sample_dict):
    component_names = ["whats_my_species", "qcquickie", "assemblatron"]
    comps = {}
    for component in component_names:
        comps[component] = datahandling.load_last_sample_component(
            str(sample_dict["_id"]), component)
    return comps

def script__test_ssi_stamper(sample, sample_yaml, sample_component):
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

    stamp_list = sample.get("stamps", [])
    stamp_list.append(stamp)
    sample["stamps"] = stamp_list
    datahandling.save_sample(sample, sample_yaml)
    return 0

script__test_ssi_stamper(snakemake.params.sample,
                        snakemake.params.sample_yaml,
                        snakemake.params.sample_component)
