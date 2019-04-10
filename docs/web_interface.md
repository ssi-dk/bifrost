# Web interface usage guide

## Description of the common pipelines

**whats_my_species:** Identifies the sample species.

**assemblatron:** Assembles the reads.

**ssi_stamper:** Checks if the sample passes SSI QC. Depends on
whats_my_species and assemblatron.

**ariba_****:** Runs mlst/*finder pipeline on the reads.
Depends on whats_my_species. *Note: ariba_virulencefinder runs only on some
species.

**min_read_check:** Checks if sample has minimum number of reads.
Not in use currently.

**qcquickie:** Runs a super fast assembly pipeline on the reads. Used to gather
data for future data analysis.

## Rerunning a component (pipeline)

(Text instructions below)

![gif reruning a component](_media/rerunning_component.gif)

- Go to the run checker
- Load your run
- Check which samples and which components for those samples you want to rerun.
- Scroll down the page to the "rerun components"
- Enter the sample and the components you want to rerun and click submit
- Everything will try to run at the same time if possible, so if you have 
  components with dependencies you should submit them after the ones they
  depend on are finished.