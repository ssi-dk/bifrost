# Adding KMA to computerome checklist

1. Once they contact us to get access to bifrost send email explaining what to do:

```

We are ready to get bifrost running on your samples. We think the easiest is for us to visit you to set up the system and discuss with you how to use it and any custom needs you might have. 

There are a few things you need before:

- A Computerome project directory to store the input and output data. To simplify the process you could have a dedicated project for bifrost but it is not required. Kim and I need access to set up the system and for debugging purposes, specially in the beginning. We are perfectly fine signing an NDA if required. Our Computerome usernames are masali and kimngku

- We need to add you as a user to our bifrost-code project so you can use the latest version of the code and the databases used by the analysis components. In order to do that, we need your Computerome username(s) so we can ask them to add you to our bifrost-code project.

- (Optional) If you have a dedicated computing node in Computerome, we can set it up so you run all your samples through it instead of submitting jobs to the full Computerome cluster. This is not required but can help in saving on computing costs if you already have a dedicated node. In that case we need what they call the Reservarion id, it is a name in the format of "system.XXXXXX" that you use when submitting using the "advres" setting. It is different from the node name and you can get it by asking Computerome tech support.

On our end, this is what we'll do:

- Create a private database for you hosted in our node, as well as a private key your bifrost installation will use. This database will include the species-specific QC parameters we use at SSI. 

- Set up a web server in our node for you to access the data you generate in the web reporter. As the development of an advanced user system is still ongoing, right now you'll have a custom web server using basic http auth with credentials that we'll send to you once we create them. 

As soon as you have the Computerome project directory and we have added you to the bifrost-code project we can head to Rigshospitalet and go through the setup and discuss the usage with you.

```

2. Add their users to ssi_10003 project

3. Create mongodb database and user account

db.createUser({ "user": "bifrost_KMA", "pwd": "random", "roles": [{role:"readWrite", db: "bifrost_KMA"}], mechanisms:["SCRAM-SHA-1"]})

add setup_date index (check which others to create)
mongo
use admin
db.auth()
use bifrost_KMA
db.sample_components.createIndex({"setup_date": 1})

- Set up web server with authentication connected to their database.


Setting up bifrost

setting up the structure

create condarc in user home
envs_dirs:
  - /home/projects/ssi_10003/apps/anaconda/envs

pkgs_dirs:
  - /home/projects/ssi_10003/apps/anaconda/pkgs


Running bifrost

cp src
ln -s samples
module load tools
module load anaconda3/4.4.0
. activate bifrost


/components/whats_my_species.yaml fix lines:
kraken_database: /home/projects/ssi_10003/data/DB/kraken/minikraken_20171019_8GB/
kraken_kmer_dist: /home/projects/ssi_10003/data/DB/kraken/minikraken_20171019_8GB/minikraken_8GB_100mers_distrib.txt

/components/analyzer.yaml
abricate_plasmidfinder_database: /home/projects/ssi_10003/data/DB/abricate/plasmidfinder_db/
abricate_resfinder_database: /home/projects/ssi_10003/data/DB/abricate/resfinder_db/
ariba_plasmidfinder_database: /home/projects/ssi_10003/data/DB/ariba/plasmidfinder/
ariba_resfinder_database: /home/projects/ssi_10003/data/DB/ariba/resfinder/


snakemake -s src/bifrost.smk --config grid=torque components=whats_my_species,analyzer,ssi_stamper group=ssi_10003
