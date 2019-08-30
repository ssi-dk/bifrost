# script for use with snakemake
import sys
import traceback
import os
from bifrostlib import datahandling


def rule__summarize_variants(input, output, sampleComponentObj, log):
    import cyvcf2
    try:
        this_function_name = sys._getframe().f_code.co_name
        sample_db, component_db = sampleComponentObj.start_rule(this_function_name, log=log)

        variants_vcf_file = str(input.variants)
        summarize_ambiguous_snp_yaml = str(output.variants_yaml)
        """
        Logic here is that you record a next_variant then compare it to the next one first. This is
        done because in the next_variant caller used (BBMap) multiple variants for the same position
        will appear on different lines. As the file is sorted positions will be beside each other
        except on the edgecase of multiple positions of 1 from new contigs which shouldn't happen.
        indels and deletions are recorded seperately and then because samples are mapped against
        themselves neither the reference or alternative represent the TRUE reference. Rather
        which ever is in a higher percent is considered the dominant one. This means that instead
        of dealing with 0% -> 100% we're dealing with 0% -> 50% as a site with 60% is really just
        saying that the base that was chosen appears to be at conflict with the mapped reads
        and since it's from the same read set and we assume the mapped are more reflective of
        reality then the alternative is instead correct at 60% and the reference in this case is
        at 40% frequency and would be recorded at 40%. These values are filled out at depth values
        of 1-100 with 100 being all depths >=100. This produces a small 100x100 data matrix which
        should be more descriptive on ambiguous sites than a value chosen at a particular cutoff.
        A value at a particular cutoff such as 10x coverage and 90% certainty could be read off
        at 10x and (100% - 90%) = 10% for the single value.
        """
        result_matrix = {}
        result_matrix["indels"] = 0
        result_matrix["deletions"] = 0
        data_matrix = {}
        for i in range(1, 101):
            data_matrix[i] = {}
            for j in range(50, -1, -1):
                data_matrix[i][j] = 0

        variant_file = cyvcf2.VCF(variants_vcf_file)
        current_variant = next(variant_file)
        current_frequency = current_variant.INFO.get("AF")
        next_frequency = 0
        for next_variant in variant_file:
            if not next_variant.is_indel and not next_variant.is_deletion:
                if current_variant.POS == next_variant.POS and current_variant.CHROM == next_variant.CHROM:
                    next_frequency = next_variant.INFO.get("AF") + current_frequency
                else:
                    next_frequency = next_variant.INFO.get("AF")
                    depth = current_variant.INFO.get("DP")
                    frequency = current_frequency
                    if depth > 100:
                        depth = 100
                    if round(frequency, 2) > 0.5:
                        frequency = 1.00 - frequency
                    frequency = int(round(frequency, 2) * 100)
                    data_matrix[depth][frequency] += 1
            elif next_variant.is_indel:
                result_matrix["indels"] += 1
            elif next_variant.is_deletion:
                result_matrix["deletions"] += 1
            current_variant = next_variant
            current_frequency = next_frequency
        # should handle last entry
        if not next_variant.is_indel and not next_variant.is_deletion:
            depth = next_variant.INFO.get("DP")
            frequency = current_frequency
            if depth > 100:
                depth = 100
            if round(frequency, 2) > 0.5:
                frequency = 1.00 - frequency
            frequency = int(round(frequency, 2) * 100)
            data_matrix[depth][frequency] += 1
        elif next_variant.is_indel:
            result_matrix["indels"] += 1
        elif next_variant.is_deletion:
            result_matrix["deletions"] += 1
    except StopIteration:
            result_matrix["error"] = "No variants"

            # back propogate the values
        variant_table = [[0 for x in range(51)] for y in range(100)]
        column_add = [0] * 51
        for i in range(100, 0, -1):
            row_add = 0
            for j in range(50, -1, -1):
                variant_table[i - 1][j] += column_add[j]
                row_add += data_matrix[i][j]
                variant_table[i - 1][j] += row_add
                column_add[j] = variant_table[i - 1][j]

        result_matrix["variant_table"] = variant_table
        datahandling.save_yaml(result_matrix, summarize_ambiguous_snp_yaml)

        sampleComponentObj.end_rule(this_function_name, log=log)
    except Exception:
        sampleComponentObj.write_log_err(log, str(traceback.format_exc()))



rule__summarize_variants(
    snakemake.input,
    snakemake.output,
    snakemake.params.sampleComponentObj,
    snakemake.log)
