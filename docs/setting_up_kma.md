# Adding KMA to computerome checklist

## 1. Once they contact us to get access to bifrost send email explaining what to do

> We are ready to get bifrost running on your samples. We think the easiest is for us to visit you to set up the system and discuss with you how to use it and any custom needs you might have.
> There are a few things you need before:
>
> - A Computerome project directory to store the input and output data. To simplify the process you could have a dedicated project for bifrost but it is not required. Kim and I need access to set up the system and for debugging purposes, specially in the beginning. We are perfectly fine signing an NDA if required. Our Computerome usernames are masali and kimngku
> - We need to add you as a user to our bifrost-code project so you can use the latest version of the code and the databases used by the analysis components. In order to do that, we need your Computerome username(s) so we can ask them to add you to our bifrost-code project.
> - (Optional) If you have a dedicated computing node in Computerome, we can set it up so you run all your samples through it instead of submitting jobs to the full Computerome cluster. This is not required but can help in saving on computing costs if you already have a dedicated node. In that case we need what they call the Reservarion id, it is a name in the format of "system.XXXXXX" that you use when submitting using the "advres" setting. It is different from the node name and you can get it by asking Computerome tech support.
>
> On our end, this is what we'll do:
>
> - Create a private database for you hosted in our node, as well as a private key your bifrost installation will use. This database will include the species-specific QC parameters we use at SSI.
> - Set up a web server in our node for you to access the data you generate in the web reporter. As the development of an advanced user system is still ongoing, right now you'll have a custom web server using basic http auth with credentials that we'll send to you once we create them.
>
> As soon as you have the Computerome project directory and we have added you to the bifrost-code project we can head to Rigshospitalet and go through the setup and discuss the usage with you.

## 2. Add their users to ssi_10003 project

## 3. MongoDB setup

### 3.1. Create mongodb database and user account

```javascript
db.createUser({ "user": "bifrost_KMA", "pwd": "random", "roles": [{role:"readWrite", db: "bifrost_KMA"}], mechanisms:["SCRAM-SHA-1"]})
```

### 3.2. Add setup_date index (check which others to create)

```bash
mongo
use admin
db.auth()
use bifrost_KMA
db.sample_components.createIndex({"setup_date": 1})
```

### 3.3. Add mongodb URI to config/keys.txt

## 4. Set up web server with authentication connected to their database

Webserver files should be in `/home/projects/KMA_project/data/bifrost/webserver/reporter/reporter.py` and run in the tmux server in our Computerome machine.

## 5. Setting up bifrost

### 5.1 Setting up the structure

```bash
mkdir -p /home/projects/KMA_project/data/bifrost/config
mkdir -p /home/projects/KMA_project/data/bifrost/input
mkdir -p /home/projects/KMA_project/data/bifrost/output
```

### 5.2 Create `.condarc` in user home (running user should do this)

```bash
envs_dirs:
  - /home/projects/ssi_10003/apps/anaconda/envs

pkgs_dirs:
  - /home/projects/ssi_10003/apps/anaconda/pkgs
```

## 6. Config

### 6.1 KMA config

Create a config.yaml file in the config directory with the following:

```yaml
grid: "torque"
group: "<groupname>"
sample_sheet: "sample_sheet.tsv"
components: "whats_my_species,analyzer,assemblatron,ssi_stamper,ariba_resfinder,ariba_mlst,ariba_plasmidfinder,ariba_virulencefinder,qcquickie,min_read_check" # available options are in components directory
type: "routine"
```

## 6. Running bifrost

See KMA_README.md for this.
