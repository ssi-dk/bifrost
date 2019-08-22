import math

def assemblatron_diff(res):
    if "assemblatron.bin_length_at_1x" in res and "assemblatron.bin_length_at_10x" in res:
        res["assemblatron.bin_length_1x_10x_diff"] = res["assemblatron.bin_length_at_1x"] - \
            res["assemblatron.bin_length_at_10x"]
    else:
        res["assemblatron.bin_length_1x_10x_diff"] = math.nan
    return res

def assemblatron_contig_diff(res):
    if "assemblatron.bin_contigs_at_1x" in res and "assemblatron.bin_contigs_at_10x" in res:
        res["assemblatron.bin_contigs_1x_10x_diff"] = res["assemblatron.bin_contigs_at_1x"] - \
            res["assemblatron.bin_contigs_at_10x"]
    else:
        res["assemblatron.bin_contigs_1x_10x_diff"] = math.nan
    return res

FUNCS = [assemblatron_diff, assemblatron_contig_diff]


COLUMNS = [
    {
        "name": "Name",
        "id": "name"
    },
    {
        "name": "Run",
        "id": "sample_sheet.run_name"
    },
    {
        "name": "Supplying_lab_feedback",
        "id": "stamps.supplying_lab_check.value",
    },
    {
        "name": "QC_action",
        "id": "stamps.ssi_stamper.value",
        # "id": "ssi_stamper.assemblatron:action"
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
        "name": "Genome_size_1x",
        "id": "properties.denovo_assembly.bin_length_at_1x"
    },
    {
        "name": "Genome_size_10x",
        "id": "properties.denovo_assembly.bin_length_at_10x"
    },
    {
        "name": "G_size_diff_1x_10",
        "id": "properties.stamper.test__denovo_assembly__genome_size_difference_1x_10x_text"
    },
    {
        "name": "Avg_coverage",
        "id": "properties.denovo_assembly.bin_coverage_at_1x"
    },
    {
        "name": "Num_contigs",
        "id": "properties.denovo_assembly.bin_contigs_at_1x"
    },
    {
        "name": "Ambiguous_sites",
        "id": "properties.denovo_assembly.snp_filter_10x_10%"
    },
    {
        "name": "Num_reads",
        "id": "properties.denovo_assembly.filtered_reads_num"
    },
    # {
    #     "name": "mlst_type",
    #     "id": "properties.mlst.strain"
    # },
    {
        "name": "Main_sp_plus_uncl",
        "id": "properties.stamper.test__species_detection__main_species_level_text"
    },
    {
        "name": "Unclassified_reads",
        "id": "properties.species_detection.percent_unclassified"
    },
    {
        "name": "DB_ID",
        "id": "_id"
    },
    {
        "name": "Failed_tests",
        "id": "ssi_stamper_failed_tests"
    }
]

plot_values = [
    {
        "name": "Genome_size_1x",
        "id": "properties.denovo_assembly.bin_length_at_1x",
        "limits": [1500000, 6000000]
    },
    {
        "name": "Genome_size_10x",
        "id": "properties.denovo_assembly.bin_length_at_10x",
        "limits": [1500000, 6000000],
        "xaxis": "x"
    },
    {
        "name": "G_size_difference_1x_10",
        "id": "properties.stamper.test__denovo_assembly__genome_size_difference_1x_10x_text",
        "limits": [0, 260000]
    },
    {
        "name": "Avg_coverage",
        "id": "properties.denovo_assembly.bin_coverage_at_1x",
        "limits": [0, 200]
    },
    {
        "name": "Contig_num_1x",
        "id": "properties.denovo_assembly.bin_contigs_at_1x",
        "limits": [0, 700]
    },
    {
        "name": "Num_reads",
        "id": "properties.denovo_assembly.filtered_reads_num",
        "limits": [1000, 8000000]
    },
    {
        "name": "Main_sp_plus_unclassified",
        "id": "properties.stamper.test__species_detection__main_species_level_text",
        "limits": [0.75, 1]
    },
    {
        "name": "Unclassified_reads",
        "id": "properties.species_detection.percent_unclassified",
        "limits": [0, 0.25]
    }
]

finder_columns = [
    {
        "name": "GENE",
        "id": "GENE"
    },
    {
        "name": "%COVERAGE",
        "id": "%COVERAGE"
    },
    {
        "name": "%IDENTITY",
        "id": "%IDENTITY"
    },
    {
        "name": "SEQUENCE",
        "id": "SEQUENCE"
    },
    {
        "name": "START",
        "id": "START"
    },
    {
        "name": "END",
        "id": "END"
    },
    {
        "name": "DATABASE",
        "id": "DATABASE"
    },
    {
        "name": "COVERAGE",
        "id": "COVERAGE"
    },
    {
        "name": "ACCESSION",
        "id": "ACCESSION"
    }
]

ROUND_COLUMNS = ["whats_my_species.percent_unclassified",
                 "assemblatron.bin_coverage_at_1x"]
