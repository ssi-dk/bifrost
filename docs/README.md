# Welcome to bifrost's documentation

Bifrost is a sharing platform/LIMS for data analysis currently focused on bacterial quality control (QC).

Bifrost is being developed by Statens Serum Institut (SSI) for tracking collections of samples at the sample level. By tracking what analyses have been run against each sample in a centralized database (DB) we can perform rapid routine analysis against all values that enter the system.

To ensure that the analyses are accurate it is important to run and track each sample through QC which is appropriate for the downstream analysis you want to perform and is thus a central component of bifrost.

A centralized DB also gives greater power in analyzing your sample collection through either direct DB interaction or via bifrost's web based graphical user interface (GUI).

Bifrost provides a platform allowing for real-time information on how a new sample compares against all existing DB samples enabling more data-based approaches for analysis.

Bifrost has been implemented at SSI and we are as of this writing in the process to rolling this platform out to all clinical microbiology labs in Denmark and will see continued development to fulfill the needs of software as infrastructure.

# Background

Bifrost is being developed with the concept that data analysis should be tracked downstream through the use of components. The idea is that summarized values are often used for downstream analysis and that we should encourage running standardized components on samples or collections of samples. This serves the purpose of both centralizing the results into a common database and encouraging the use of standardized tools for analysis, while still allowing the flexibility of using personalized tools. The first components for bifrost have all been developed based on quality control of bacterial WGS samples. If the resulting data is standardized then visual analysis tools can be built on top of the summarized data and pool information from all samples that have also run the same component and additional metrics can be discovered more easilly off the pooled datasets.

Bifrost originated from Statens Serum Institut (SSI) to handle the increase of whole genome sequence (WGS) samples of bacterial isolates that are processed year by year. Currently over 10,000 samples of varying species are processed each year. The increase in samples for outbreak, surveillance and project use required a common QC platform for the re-use of samples in downstream analysis. Over time the common thread of using sequences over and over became more apparent. Tracking project to project or within group specialities would work for a while but would not maximize the value of the data by being able to aggregate across all of our data. 

On top of this greater data sharing will occur between the microbiology labs in Denmark and solutions on how to facilitate this transfer was also needed. The reality of the situation is that there is plenty of informtion from different locations that just cannot be shared fully, at least not currently.
