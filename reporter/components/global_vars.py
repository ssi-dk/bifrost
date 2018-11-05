import collections
import math

PLOTS = collections.OrderedDict()

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


FUNCS = [assemblatron_diff, assemblatron_contig_diff]
