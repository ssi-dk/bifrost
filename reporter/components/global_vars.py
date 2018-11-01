import collections
import math

PLOTS = collections.OrderedDict()

# PLOTS["qcquickie.bin_length_at_1x"] = {
#     "projection": "qcquickie.bin_length_at_1x"
# }
# PLOTS["qcquickie.bin_length_at_10x"] = {
#     "projection": "qcquickie.bin_length_at_10x"
# }
# PLOTS["qcquickie.bin_length_at_25x"] = {
#     "projection": "qcquickie.bin_length_at_25x"
# }

def qcquickie_diff(res):
    if "qcquickie.bin_length_at_1x" in res and "qcquickie.bin_length_at_25x" in res:
        res["qcquickie.bin_length_1x_25x_diff"] = res["qcquickie.bin_length_at_1x"] - \
            res["qcquickie.bin_length_at_25x"]
    return res


# PLOTS["qcquickie.bin_coverage_at_1x"] = {
#     "projection": "qcquickie.bin_coverage_at_1x"
# }
# PLOTS["qcquickie.bin_coverage_at_10x"] = {
#     "projection": "qcquickie.bin_coverage_at_10x"
# }
# PLOTS["qcquickie.bin_coverage_at_25x"] = {
#     "projection": "qcquickie.bin_coverage_at_25x"
# }
# PLOTS["qcquickie.bin_length_1x_25x_diff"] = {
#     "projection": "qcquickie.bin_length_1x_25x_diff",
#         "func": qcquickie_diff
#     }
# PLOTS["qcquickie.bin_contigs_at_1x"] = {
#     "projection": "qcquickie.bin_contigs_at_1x"
# }
# PLOTS["qcquickie.bin_contigs_at_10x"] = {
#     "projection": "qcquickie.bin_contigs_at_10x"
# }
# PLOTS["qcquickie.bin_contigs_at_25x"] = {
#     "projection": "qcquickie.bin_contigs_at_25x"
# }
# PLOTS["qcquickie.N50"] = {
#     "projection": "qcquickie.N50"
# }
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
    if "assemblatron.bin_length_at_1x" in res and "assemblatron.bin_length_at_10x" in res:
        res["assemblatron.bin_length_1x_10x_diff"] = res["assemblatron.bin_length_at_1x"] - \
            res["assemblatron.bin_length_at_10x"]
    else:
        res["assemblatron.bin_length_1x_10x_diff"] = math.nan
    return res

PLOTS["assemblatron.bin_length_1x_10x_diff"] = {
    "projection": "assemblatron.bin_length_1x_10x_diff",
        "func": assemblatron_diff
    }

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

def assemblatron_contig_diff(res):
    if "assemblatron.bin_contigs_at_1x" in res and "assemblatron.bin_contigs_at_10x" in res:
        res["assemblatron.bin_contigs_1x_10x_diff"] = res["assemblatron.bin_contigs_at_1x"] - \
            res["assemblatron.bin_contigs_at_10x"]
    else:
        res["assemblatron.bin_contigs_1x_10x_diff"] = math.nan
    return res


PLOTS["assemblatron.bin_contigs_1x_10x_diff"] = {
    "projection": "assemblatron.bin_contigs_1x_10x_diff",
    "func": assemblatron_contig_diff
}

PLOTS["assemblatron.bin_contigs_at_25x"] = {
    "projection": "assemblatron.bin_contigs_at_25x"
}
PLOTS["assemblatron.N50"] = {
    "projection": "assembly.N50"
}

DEFAULT_PLOT = "assemblatron.bin_length_at_1x"
