# Development

Bifrost is meant to be an ambitious project for data analysis tracking. To accomplish this the following development projects are being explored in future development:

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
Since we code in python and have been working with snakemake there has been no local interest to develop handling for other workflows. In fact there has been more talk in doing the opposite and adjusting our launcher to also be fully integrated in snakemake

## More automated testing (internal)

### Idea
With the size of the system scope the developers have adopted a devops approach for managing it. To adequately handle this more test cases need to be created and continually tested against the code base to free up time from testing development changes. While this is specifically applied to our internal platform we are poositive many of the test cases can be pushed towards the community at large.

### Technical
So this is already in progress and part of it is handling it within our Azure DevOps framework to maximize the utilization of that. Part of it is to still learn more about the Azure DevOps framework and how it's supposed to work properly and for our development team. The goal is everytime we have a bug a test gets built up to ensure it gets auto tested to ensure we know about it when it occurs again. More bugs equals more tests and conceptually this will decrease as we develop more for this. Because of the nature of our platform we have to ensure we can run those tests appropriately and which builds we run it on.

Better handling of metadata
Examination of existing data to better facilitate background information
Addition of Tags for sample/run
Private/Public databases
Autorun of outdated components
Refactoring of bifrost launcher
API development
Import/Export of DB data
Run components
Anonymization of data
Workspaces/Projects (on run)
