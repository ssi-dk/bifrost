import collections

columns = ['name', 'supplying_lab', 'run_name', '_id',
            'input_read_status', 'emails', 'user',
            'R1_location', 'R2_location',
            'provided_species',
            'qcquickie_name_classified_species_1',
            'qcquickie_percent_classified_species_1',
            'qcquickie_name_classified_species_2',
            'qcquickie_percent_classified_species_2',
            'qcquickie_percent_unclassified',
            'qcquickie_bin_length_at_1x',
            'qcquickie_bin_length_at_10x',
            'qcquickie_bin_length_at_25x',
            'qcquickie_bin_length_1x_25x_diff',
            'qcquickie_bin_coverage_at_1x',
            'qcquickie_bin_coverage_at_10x',
            'qcquickie_bin_coverage_at_25x',
            'qcquickie_bin_contigs_at_1x',
            'qcquickie_bin_contigs_at_10x',
            'qcquickie_N50', 'qcquickie_N75',
            'qcquickie_snp_filter_deletions',
            'qcquickie_snp_filter_indels',
            'qcquickie_snp_filter_10x_10%',
            'assembly_bin_length_at_1x',
            'assembly_bin_length_at_10x',
            'assembly_bin_length_at_25x',
            'assembly_bin_length_1x_25x_diff',
            'assembly_bin_coverage_at_1x',
            'assembly_bin_coverage_at_10x',
            'assembly_bin_coverage_at_25x',
            'assembly_bin_contigs_at_1x',
            'assembly_bin_contigs_at_10x',
            'assembly_N50', 'assembly_N75',
            'assembly_snp_filter_deletions',
            'assembly_snp_filter_indels',
            'assembly_snp_filter_10x_10%', 'comments']


PLOTS = collections.OrderedDict()

PLOTS["qcquickie_bin_length_at_1x"] = {
        "projection": "qcquickie.summary.bin_length_at_1x"
    }
PLOTS["qcquickie_bin_length_at_10x"] = {
        "projection": "qcquickie.summary.bin_length_at_10x"
    }
PLOTS["qcquickie_bin_length_at_25x"] = {
        "projection": "qcquickie.summary.bin_length_at_25x"
    }
PLOTS["qcquickie_bin_length_1x_25x_diff"] = {
        "projection": "value",
        "aggregate": [
            {
                "$project": {
                    "sample.name": 1,
                    "qcquickie.summary.name_classified_species_1": 1,
                    "value": {
                        "$subtract": [
                            "$qcquickie.summary.bin_length_at_1x",
                            "$qcquickie.summary.bin_length_at_25x"
                        ]
                    }
                },
            },
        ]
    }
PLOTS["qcquickie_bin_contigs_at_1x"] = {
        "projection": "qcquickie.summary.bin_contigs_at_1x"
    }
PLOTS["qcquickie_bin_contigs_at_10x"] = {
        "projection": "qcquickie.summary.bin_contigs_at_10x"
    }
PLOTS["qcquickie_bin_contigs_at_25x"] = {
        "projection": "qcquickie.summary.bin_contigs_at_25x"
    }
PLOTS["qcquickie_N50"] = {
        "projection": "qcquickie.quast/report_tsv/N50"
    }
PLOTS["qcquickie_N75"] = {
        "projection": "qcquickie.quast/report_tsv/N75"
    }
PLOTS["assembly_bin_length_at_1x"] = {
        "projection": "assembly.summary.bin_length_at_1x"
    }
PLOTS["assembly_bin_length_at_10x"] = {
        "projection": "assembly.summary.bin_length_at_10x"
    }
PLOTS["assembly_bin_length_at_25x"] = {
        "projection": "assembly.summary.bin_length_at_25x"
    }
PLOTS["assembly_bin_length_1x_25x_diff"] = {
        "projection": "value",
        "aggregate": [
            {
                "$project": {
                    "sample.name": 1,
                    "qcquickie.summary.name_classified_species_1": 1,
                    "value": {
                        "$subtract": [
                            "$assembly.summary.bin_length_at_1x",
                            "$assembly.summary.bin_length_at_25x"
                        ]
                    }
                },
            },
        ]
    }
PLOTS["assembly_bin_contigs_at_1x"] = {
        "projection": "assembly.summary.bin_contigs_at_1x"
    }
PLOTS["assembly_bin_contigs_at_10x"] = {
        "projection": "assembly.summary.bin_contigs_at_10x"
    }
PLOTS["assembly_bin_contigs_at_25x"] = {
        "projection": "assembly.summary.bin_contigs_at_25x"
    }
PLOTS["assembly_N50"] = {
        "projection": "assembly.quast/report_tsv/N50"
    }
PLOTS["assembly_N75"] = {
        "projection": "assembly.quast/report_tsv/N75"
    }

DEFAULT_PLOT = 'qcquickie_bin_length_at_1x'
