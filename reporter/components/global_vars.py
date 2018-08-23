import collections
import math

columns = ["name", "supplying_lab", "run_name", "_id",
            "input_read_status", "emails", "user",
            "R1_location", "R2_location",
            "provided_species",
            "qcquickie.name_classified_species_1",
            "qcquickie.percent_classified_species_1",
            "qcquickie.name_classified_species_2",
            "qcquickie.percent_classified_species_2",
            "qcquickie.percent_unclassified",
            "qcquickie.bin_length_at_1x",
            "qcquickie.bin_length_at_10x",
            "qcquickie.bin_length_at_25x",
            "qcquickie.bin_length_1x_25x_diff",
            "qcquickie.bin_coverage_at_1x",
            "qcquickie.bin_coverage_at_10x",
            "qcquickie.bin_coverage_at_25x",
            "qcquickie.bin_contigs_at_1x",
            "qcquickie.bin_contigs_at_10x",
            "qcquickie.N50",
            "qcquickie.snp_filter_deletions",
            "qcquickie.snp_filter_indels",
            "qcquickie.snp_filter_10x_10%",
            "assemblatron.bin_length_at_1x",
            "assemblatron.bin_length_at_10x",
            "assemblatron.bin_length_at_25x",
            "assemblatron.bin_length_1x_25x_diff",
            "assemblatron.bin_coverage_at_1x",
            "assemblatron.bin_coverage_at_10x",
            "assemblatron.bin_coverage_at_25x",
            "assemblatron.bin_contigs_at_1x",
            "assemblatron.bin_contigs_at_10x",
            "assemblatron.N50",
            "assemblatron.snp_filter_deletions",
            "assemblatron.snp_filter_indels",
            "assemblatron.snp_filter_10x_10%", "comments"]


PLOTS = collections.OrderedDict()

PLOTS["qcquickie.bin_length_at_1x"] = {
    "projection": "qcquickie.bin_length_at_1x"
}
PLOTS["qcquickie.bin_length_at_10x"] = {
    "projection": "qcquickie.bin_length_at_10x"
}
PLOTS["qcquickie.bin_length_at_25x"] = {
    "projection": "qcquickie.bin_length_at_25x"
}

def qcquickie_diff(res):
    if "qcquickie.bin_length_at_1x" in res and "qcquickie.bin_length_at_25x" in res:
        res["qcquickie.bin_length_1x_25x_diff"] = res["qcquickie.bin_length_at_1x"] - \
            res["qcquickie.bin_length_at_25x"]
    return res


PLOTS["qcquickie.bin_coverage_at_1x"] = {
    "projection": "qcquickie.bin_coverage_at_1x"
}
PLOTS["qcquickie.bin_coverage_at_10x"] = {
    "projection": "qcquickie.bin_coverage_at_10x"
}
PLOTS["qcquickie.bin_coverage_at_25x"] = {
    "projection": "qcquickie.bin_coverage_at_25x"
}
PLOTS["qcquickie.bin_length_1x_25x_diff"] = {
    "projection": "qcquickie.bin_length_1x_25x_diff",
        "func": qcquickie_diff
    }
PLOTS["qcquickie.bin_contigs_at_1x"] = {
    "projection": "qcquickie.bin_contigs_at_1x"
}
PLOTS["qcquickie.bin_contigs_at_10x"] = {
    "projection": "qcquickie.bin_contigs_at_10x"
}
PLOTS["qcquickie.bin_contigs_at_25x"] = {
    "projection": "qcquickie.bin_contigs_at_25x"
}
PLOTS["qcquickie.N50"] = {
    "projection": "qcquickie.N50"
}
PLOTS["assemblatron.bin_length_at_1x"] = {
    "projection": "assemblatron.bin_length_at_1x"
}
PLOTS["assemblatron.bin_length_at_10x"] = {
    "projection": "assemblatron.bin_length_at_10x"
}
PLOTS["assemblatron.bin_length_at_25x"] = {
    "projection": "assemblatron.bin_length_at_25x"
}
def assemblatron_diff(res):
    if "assemblatron.bin_length_at_1x" in res and "assemblatron.bin_length_at_25x" in res:
        res["assemblatron.bin_length_1x_25x_diff"] = res["assemblatron.bin_length_at_1x"] - \
            res["assemblatron.bin_length_at_25x"]
    else:
        res["assemblatron.bin_length_1x_25x_diff"] = math.nan
    return res


PLOTS["assemblatron.bin_coverage_at_1x"] = {
    "projection": "assemblatron.bin_coverage_at_1x"
}
PLOTS["assemblatron.bin_coverage_at_10x"] = {
    "projection": "assemblatron.bin_coverage_at_10x"
}
PLOTS["assemblatron.bin_coverage_at_25x"] = {
    "projection": "assemblatron.bin_coverage_at_25x"
}
PLOTS["assemblatron.bin_contigs_at_1x"] = {
    "projection": "assemblatron.bin_contigs_at_1x"
}
PLOTS["assemblatron.bin_contigs_at_10x"] = {
    "projection": "assemblatron.bin_contigs_at_10x"
}
PLOTS["assemblatron.bin_contigs_at_25x"] = {
    "projection": "assemblatron.bin_contigs_at_25x"
}
PLOTS["assemblatron.N50"] = {
    "projection": "assembly.N50"
}

DEFAULT_PLOT = "qcquickie.bin_length_at_1x"
