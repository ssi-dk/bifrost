# Development

Bifrost is meant to be an ambitious project for data analysis tracking. To accomplish this the following development projects are being explored in future development:


## More automated testing (internal)

### Idea
With the size of the system scope the developers have adopted a devops approach for managing it. To adequately handle this more test cases need to be created and continually tested against the code base to free up time from testing development changes. While this is specifically applied to our internal platform we are poositive many of the test cases can be pushed towards the community at large.

### Technical
So this is already in progress and part of it is handling it within our Azure DevOps framework to maximize the utilization of that. Part of it is to still learn more about the Azure DevOps framework and how it's supposed to work properly and for our development team. The goal is everytime we have a bug a test gets built up to ensure it gets auto tested to ensure we know about it when it occurs again. More bugs equals more tests and conceptually this will decrease as we develop more for this. Because of the nature of our platform we have to ensure we can run those tests appropriately and which builds we run it on.

## Examination of existing data to better facilitate background information

### Idea
With our institute alone we have tens of thousands of samples already ran through the same set of components. Unsurprisingly new insights on better metrics can be made from this existing data set especially when we store all that information into a centralized DB for many values. The idea is to adjust the settings we use for new QC metrics based on the data we already have. This can be shown with multiple distributions showing the pass distribution overlaid with the fail distribution to help better educate the decision going forward. Also you can figure out how your samples compare to various optimal and make adjustments in the lab to optimize either in quality or cost. Ideally this knowledge can be exported to be available to others so that someone who is sequencing their first sample can then use the metrics we provide to see how their sample compares which is of course very nice for people lacking that information or comparing lab workflows between institutions.

### Technical
This is a pretty broad idea and at first will be reviewing the analytics of our samples to see if our metrics are being well used. Additional learning tools can be used to also see if metrics being recorded but not used for QC are better at determining how well the sample will turn out. Outside of this with the information we can go a step forward and use the samples which are good to optimize aspects such as databases we use for checking species. 


## Component compartmentalization (Continiuum/Docker)

### Idea
Right now the required programs for each component of bifrost is installed against the base bifrost environment. This unfortunately can cause collisions between different pipelines and makes managing the seperation of components more difficult in order to ensure that each component has the right version of database and programs 
The plan is to have it that each component is launched through it's own container which contains all programs and resources already loaded into it and be provided a specific version for each container. This will both simplify the management of components as they'll be self contained, and mean that sharing a container will mean the same programs, resources, and environment for whoever uses it. 

### Technical
In order to handle compartmentalizaation the plan is to create a dockerimg for each component. Each component will have all the hooks required to function on an individual samples without any other information so that running a component and passing the inputs should yield the expected output. As there's a documented input/output for the database in each component all checks should be able to be made prior to running and output. Snakemake already offers the use of continuum images for it though some settings have to be set to load these from a centralized source instead of having them created on the fly. By including databases within the image we can really ensure the sameness of pacakges for accreditation purposes. This does however mean that new images are created routinely according to the update in the background databases. Some thoughts on how frequent to do this also need to be considered as well as the naming convention for this.


## New bacterial WGS component development

### Idea
New components are being developed to facilitate more analysis. Conceptually this goes many ways, most directly it's development of new components that further extend development related to bacterial WGS, for example running cgMLST on each sample and additional finders. Currently I'm switching the ariba pipelines to utilize kma for internal use and issues with the queue system that ariba has. 

### Technical
With how bifrost is set up anyone can make a component and run it against their sample. What we want to avoid is multiple similar components. In developing new components we want to provide users with flexibility to match their workflows but subtly encourage everyone to the same pipeline set. The hope on avoiding this is through the background data that components provide. If we as an institute only use a subset of components and provide the large number of samples as background for others using the component then that can have meaningful value. For example seeing how the QC value of your new sample compares when there is an existing background of 1000's of samples from our institute instead of seeing it on it's own provides a degree of context.


## Component development for non bacterial WGS

### Idea
On a more broad perspective components can be made for non bacterial WGS. For example microbiome samples, meta data samples, Human samples, etc. Because the concept of bifrost is to track analysis it can of course be applied to other analysis types, it's just that we started with bacterial WGS QC because it matched our usecase. 

### Technical
Converting other pipelines into bifrost was designed to have minimal new items when compared to snakemake. Anything that your institute runs against common pipelines. Therefore, the growth of components also has to be put towards whoever will utilize them. 


## Allow other workflows (nextflow) to also work with bifrost or go full snakemake

### Idea
Right now the current bifrost.smk is merely a launcher of snakemake jobs. It can be shifted to work with other workflows to make it more flexible, the new workflows will have to integrate with the APIs in order to communicate to the DB appropriately. Alternatively embracing snakemake fully can be beneficial for testing by adjusting the launcher to also use snakemake and integrating more standardly into snakemake structure

### Technical
Since we code in python and have been working with snakemake there has been no local interest to develop handling for other workflows. In fact there has been more talk in doing the opposite and adjusting our launcher to also be fully integrated in snakemake. Since there's currently no development outside our own I'd lean towards more snakemake integration to better utilize there tools for such things as queue handling and resource management. It would also be possible to wrap other commands in snakemake which is less ideal but not too hard. This will likely be a bit of wait and see as it's lower in development priority.


## Run components

### Idea
Currently all components are implemented on a per sample basis. The idea of run components goes to the idea that a run is really a group of samples and that many analysis only have contextual meaning as a grouping. For example in bacterial WGS doing a phylogenetic tree is dependant on the samples included. This would be an example of a run component. 

### Technical
The database pretty much already supports this but the components have to be made and the interfaces for them as well as well as the launcher for handling them. Ideally there is also integration between run components and some sample components within the run component. Once this gains more momentum we'll review the DB naming (as runs is less descriptive) and adjust the launcher to handle this. I think creating the components themselves will not be too challenging (assuming we're converting existing workflows like SNP based clustering)


## Better handling of metadata

### Idea
Currently the analysis are all based off pipelines that have been executed against the input or can also be considered genetic information. We have information on lots of samples in terms of meta data. Also at the core concept of a sample a sample isn't just a pair of reads. As long as a component can be created for it (either run or sample) we can also have samples which are nothing but epi data to provide context (ie for Surveillance purpose when only a subset of samples are sequenced and epi data exists for other likely samples)

### Technical
While it's currently possible to store this information into the sample itself a more comprehensive approach will likely need to be taken when more epi information gets included. Also because not all epi information applies to individual samples much like runs. This will evolve when we start a project which will address this and there is a likely possibility a new DB element for epidata is created. 

## Addition of Tags for sample/run

### Idea
Tags are essemtially a piece of meta data that we could assign to a sample or run to help organize information. It's a simple binary piece of metadata that would allow for easier quering of relative information.

### Technical
Implementing tags is mostly about implementing it into the interface. Organizing it for now is simply adding a value in the DB for tags which stores a list of applicable tabs that could be searched on. Then users could tag their items (or programs such as it's resistance profile) and the tags could be processed for better querying.

## Private/Public databases

### Idea
The idea here originated with sharing with other institutions on our DB structure. Because of how it works conceptually we could think we have 2 DB's at each institute a private and a public one. The private could contain more details than the public one or finer detail (e.g. gps coords vs country coords). The public one could then be shared within a trusted circle or freely while the private one would remain for personal use.

### Technical
So the idea is because we are using mongoDB on our implementation it offers options for sharing data via a sharded system. Information could be shared without actually having to put all information in one location (as long as you could connect to each source). You'd also have a mapping between the public and private so if inquiries are made on partial data additional information could be provided if desired. I don't think two seperated DB are actually required but it might more sense to segment it that way for conceptualization and being sure you put data in public or private.


## Autorun of outdated components

### Idea
Currently we're on the same set of components however the expectation is a newer better software will eventually come and that our component will need an update. When this happens some analysis will only work with all samples having been run on the same version of component. In cases like this it would be nice to have a background process which automatically ran the backlog of samples to upgrade them.

### Technical
This is dependant on the environment but currently with ours we can monitor our queue sizes for free space on our system and run jobs in the background based on queries. Once we develop a new component that requires this backlog we will likely implement something for our solution and hopefully that proves useful for others if for nothing else the queries and monitoring of components. This may also have to occur with notable database updates. Some thought is being put in to handle databases which are continually updated.


## Refactoring of bifrost launcher

### Idea
The launcher for bifrost doesn't utilize all the tools of snakemake very well and while functional for our implementations could have issues with other environments. It also just needs a general refactoring with the increase in knowledge we have now for addressing code.

### Technical
The current plan is to incorporate checkpoints into components and remake the bifrost launcher in snakemake incorporating the new components so that everything can run 100% in snakemake which helps us leverage the benefits of snakemake more so. Other changes for the launcher would revolve around being able to incorporate run components into it seemlessly and monitoring progress of steps.


## API development

### Idea
Right now we're making calls to datahandling which are really working directly against the yaml files and not the database. This should be expanded to a proper API which can handle DB requests consistently and without direct interaction.

### Technical
The plan is to develop a full rest api for handling data interactions which should also shore up how the system is accessed more securely. It should also open up more options for developing components and transferring data for the system. Ideally all requests would go through this.


## Import/Export of DB data

### Idea
With all the planned sharring we need processes set up to export the data from one database while doing all the proper handshakes and importing it at a different location (automatically or manually).

### Technical
This will be built with the API development so that data can be tranferred through this system. It will also likely take part of the data anonymization in this development. The idea is of course to be able to transfer data between different systems while maintaing proper hooks. If it creates new samples that's fairly easy but we also need ways to map up and update existing entries.

## Anonymization of data

### Idea
Sharing of data is not always able to be done at a full level. Sometimes names or metadata contain sensative data othertimes some information only makes sense in internal projects. The idea is we can mask some or all values of a sample out for sharing purposes so that something can be provided instead of nothing. 

### Technical
Technically we're looking at a few things. One is a table mapping the changes for each anonymization to the local ID's so that they can be back tracked in case additional information is needed in the future. Ideally when info is shared contact information is also shared for this process. This means you can share the same data different ways to different people. This may also be applicable on an internal sense to avoid naming samples and instead tracking them on their id's or an internally unique naming scheme so samples can be exported with different names if need be. Regardless the samples can be shared with masked information and imported into a new implementation to be used for other processing.

## Workspaces/Projects (on run)

### Idea
So right now everything is handled in term of raw data runs however runs are really synonymous with a group of samples. This could be a way to organize projects or workspaces where samples are being changed in and out. Interfaces for this need to be done to make it workable.

### Technical
This is similar to tags but handles more groupiings in a more common manner and flows into run components much more nicely. Having interfaces to allow this will then allow people to manage there samples and run downstream run components by knowing where everything is and how it should be organized.
