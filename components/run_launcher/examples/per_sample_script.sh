# Per sample script example
# This example will assuming paths are set up, run each sample against the 3
# different components (min_read_check v2.0.6, whats_my_species v2.0.6, assemblatron v2.0.6)
# using singularity, you can add other commands here including grid engine.
#
# NOTE_1: You can access $sample and $run variables with $run.sub_value.another_sub_value.etc
# this shouldn't have issues as mongoDB doesn't allow .'s in key names so we can
# parse on that value.
echo "Running $sample.name from $run.name";
BIFROST_RAW_DATA_MNT="/raw_data/mnt";
BIFROST_PIPELINE_TOOLS="/tools/singularity";
mkdir $sample.name;
cd $sample.name;
singularity run -B \
$BIFROST_RAW_DATA_MNT,\
$sample.properties.paired_reads.summary.data[0],\
$sample.properties.paired_reads.summary.data[1] \
$BIFROST_PIPELINE_TOOLS/bifrost-min_read_check_2.0.6.sif \
-id $sample._id;
singularity run -B \
$BIFROST_RAW_DATA_MNT,\
$sample.properties.paired_reads.summary.data[0],\
$sample.properties.paired_reads.summary.data[1] \
$BIFROST_PIPELINE_TOOLS/bifrost-whats_my_species_2.0.6.sif \
-id $sample._id;
singularity run -B \
$BIFROST_RAW_DATA_MNT,\
$sample.properties.paired_reads.summary.data[0],\
$sample.properties.paired_reads.summary.data[1] \
$BIFROST_PIPELINE_TOOLS/bifrost-assemblatron_2.0.6.sif \
-id $sample._id;
cd ..;
echo "Done $sample.name from $run.name";
