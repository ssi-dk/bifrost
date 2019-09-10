import math

def assemblatron_diff(res):
    if "properties.denovo_assembly.summary.bin_length_at_1x" in res and "properties.denovo_assembly.summary.bin_length_at_10x" in res:
        res["properties.denovo_assembly.summary.bin_length_1x_10x_diff"] = res["properties.denovo_assembly.summary.bin_length_at_1x"] - \
            res["properties.denovo_assembly.summary.bin_length_at_10x"]
    else:
        res["properties.denovo_assembly.summary.bin_length_1x_10x_diff"] = math.nan
    return res

def assemblatron_contig_diff(res):
    if "properties.denovo_assembly.summary.bin_contigs_at_1x" in res and "properties.denovo_assembly.summary.bin_contigs_at_10x" in res:
        res["properties.denovo_assembly.summary.bin_contigs_1x_10x_diff"] = res["properties.denovo_assembly.summary.bin_contigs_at_1x"] - \
            res["properties.denovo_assembly.summary.bin_contigs_at_10x"]
    else:
        res["properties.denovo_assembly.summary.bin_contigs_1x_10x_diff"] = math.nan
    return res

FUNCS = [assemblatron_diff, assemblatron_contig_diff]


COLUMNS = [
    {
        "name": "Run",
        "id": "sample_sheet.run_name"
    },
    {
        "name": "Name",
        "id": "name"
    },
    {
        "name": "Supplying_lab_feedback",
        "id": "stamps.supplying_lab_check.value",
    },
    {
        "name": "QC_action",
        "id": "properties.stamper.summary.stamp.value",
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
        "id": "properties.species_detection.summary.provided_species"
    },
    {
        "name": "Detected Species",
        "id": "properties.species_detection.summary.detected_species"
    },
    {
        "name": "Genome_size_1x",
        "id": "properties.denovo_assembly.summary.bin_length_at_1x"
    },
    {
        "name": "Genome_size_10x",
        "id": "properties.denovo_assembly.summary.bin_length_at_10x"
    },
    {
        "name": "G_size_diff_1x_10",
        "id": "properties.stamper.summary.test__denovo_assembly__genome_size_difference_1x_10x.value"
    },
    {
        "name": "Avg_coverage",
        "id": "properties.denovo_assembly.summary.bin_coverage_at_1x"
    },
    {
        "name": "Num_contigs",
        "id": "properties.denovo_assembly.summary.bin_contigs_at_1x"
    },
    {
        "name": "Ambiguous_sites",
        "id": "properties.denovo_assembly.summary.snp_filter_10x_10%"
    },
    {
        "name": "Num_reads",
        "id": "properties.denovo_assembly.summary.filtered_reads_num"
    },
    {
        "name": "mlst_type",
        "id": "properties.mlst.summary.strain"
    },
    {
        "name": "Main_sp_plus_uncl",
        "id": "properties.stamper.summary.test__species_detection__main_species_level.value"
    },
    {
        "name": "Unclassified_reads",
        "id": "properties.species_detection.summary.percent_unclassified"
    },
    {
        "name": "File path",
        "id": "path"
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
        "id": "properties.denovo_assembly.summary.bin_length_at_1x",
        "limits": [1500000, 6000000]
    },
    {
        "name": "Genome_size_10x",
        "id": "properties.denovo_assembly.summary.bin_length_at_10x",
        "limits": [1500000, 6000000],
        "xaxis": "x"
    },
    {
        "name": "G_size_difference_1x_10",
        "id": "properties.stamper.summary.test__denovo_assembly__genome_size_difference_1x_10x.value",
        "limits": [0, 260000]
    },
    {
        "name": "Avg_coverage",
        "id": "properties.denovo_assembly.summary.bin_coverage_at_1x",
        "limits": [0, 200]
    },
    {
        "name": "Contig_num_1x",
        "id": "properties.denovo_assembly.summary.bin_contigs_at_1x",
        "limits": [0, 700]
    },
    {
        "name": "Num_reads",
        "id": "properties.denovo_assembly.summary.filtered_reads_num",
        "limits": [1000, 8000000]
    },
    {
        "name": "Main_sp_plus_unclassified",
        "id": "properties.stamper.summary.test__species_detection__main_species_level.value",
        "limits": [0.75, 1]
    },
    {
        "name": "Unclassified_reads",
        "id": "properties.species_detection.summary.percent_unclassified",
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

ROUND_COLUMNS = ["properties.species_detection.summary.percent_unclassified",
                 "properties.denovo_assembly.summary.bin_coverage_at_1x"]

expected_results = ["resistance", "mlst", "plasmid",
                    "virulence"]  # expected categories in sample report
