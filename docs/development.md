# Development

Bifrost is meant to be an ambitious project for data analysis tracking. To accomplish this the following development projects are being explored in future development:

## Component compartmentalization (Continiuum/Docker)

### Idea
Right now the required programs for each component of bifrost is installed against the base bifrost environment. This unfortunately can cause collisions between different pipelines and makes managing the seperation of components more difficult in order to ensure that each component has the right version of database and programs 
The plan is to have it that each component is launched through it's own container which contains all programs and resources already loaded into it and be provided a specific version for each container. This will both simplify the management of components as they'll be self contained, and mean that sharing a container will mean the same programs, resources, and environment for whoever uses it. 

### Technical
In order to handle compartmentalizaation the plan is to create a dockerimg for each component. Each component will have all the hooks required to function on an individual samples without any other information so that running a component and passing the inputs should yield the expected output. As there's a documented input/output for the database in each component all checks should be able to be made prior to running and output. Snakemake already offers the use of continuum images for it though some settings have to be set to load these from a centralized source instead of having them created on the fly. By including databases within the image we can really ensure the sameness of pacakges for accreditation purposes. This does however mean that new images are created routinely according to the update in the background databases. Some thoughts on how frequent to do this also need to be considered as well as the naming convention for this.

## More automated testing (internal)

### Idea
With the size of the system scope the developers have adopted a devops approach for managing it. To adequately handle this more test cases need to be created and continually tested against the code base to free up time from testing development changes. While this is specifically applied to our internal platform we are poositive many of the test cases can be pushed towards the community at large.

Better handling of metadata
Examination of existing data to better facilitate background information
Addition of Tags for sample/run
Private/Public databases

Autorun of outdated components
Workspaces/Projects (on run)
Refactoring of bifrost launcher
API development
Import/Export of DB data
Run components
Component development
Anonymization of data