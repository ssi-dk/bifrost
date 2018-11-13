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


COLUMNS = [
    {
        "name": "Name",
        "id": "name"
    },
    {
        "name": "QC_action",
        "id": "testomatic.assemblatron:action"
    },
    {
        "name": "Comments",
        "id": "sample_sheet.Comments"
    },
    {
        "name": "Supplying_lab",
        "id": "sample_sheet.group"
    },
    {
        "name": "Provided_Species",
        "id": "properties.provided_species"
    },
    {
        "name": "Detected Species",
        "id": "properties.detected_species"
    },
    {
        "name": "QC_Genome_size_1x",
        "id": "testomatic.assemblatron:1xgenomesize"
    },
    {
        "name": "Genome_size_1x",
        "id": "assemblatron.bin_length_at_1x"
    },
    {
        "name": "Genome_size_10x",
        "id": "testomatic.assemblatron:10xgenomesize"
    },
    {
        "name": "QC_Genome_size_10x",
        "id": "assemblatron.bin_length_at_1x"
    },
    {
        "name": "QC_G_size_difference_1x_10",
        "id": "QC_testomatic.assemblatron:1x10xsizediff"
    },
    {
        "name": "G_size_difference_1x_10",
        "id": "testomatic.assemblatron:1x10xsizediff"
    },
    {
        "name": "QC_Avg_coverage",
        "id": "testomatic.assemblatron:avgcoverage"
    },
    {
        "name": "Avg_coverage",
        "id": "assemblatron.bin_coverage_at_1x"
    },
    {
        "name": "Num_contigs",
        "id": "assemblatron.bin_contigs_at_1x"
    },
    {
        "name": "Ambiguous_sites",
        "id": "assemblatron.snp_filter_10x_10%"
    },
    {
        "name": "QC_Num_reads",
        "id": "testomatic.assemblatron:numreads"
    },
    {
        "name": "Num_reads",
        "id": "assemblatron.filtered_reads_num"
    },
    {
        "name": "mlst",
        "id": "analyzer.mlst_report"
    },
    {
        "name": "QC_Main_species_read_percent",
        "id": "QC_testomatic.whats_my_species:minspecies"
    },
    {
        "name": "Main_species_read_percent",
        "id": "testomatic.whats_my_species:minspecies"
    },
    {
        "name": "QC_Unclassified_reads",
        "id": "testomatic.whats_my_species:maxunclassified"
    },
    {
        "name": "Unclassified_reads",
        "id": "whats_my_species.percent_unclassified"
    },
    {
        "name": "QC_Detected_species_in_DB",
        "id": "QC_testomatic.whats_my_species:nosubmitted"
    },
    {
        "name": "Detected_species_in_DB",
        "id": "testomatic.whats_my_species:nosubmitted"
    },
    {
        "name": "QC_Submitted_sp_is_same_as_detected",
        "id": "QC_testomatic.whats_my_species:submitted==detected"
    },
    {
        "name": "Submitted_sp_is_same_as_detected",
        "id": "testomatic.whats_my_species:submitted==detected"
    },
    {
        "name": "DB_ID",
        "id": "_id"
    },

]
