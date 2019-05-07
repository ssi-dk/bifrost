# Information
# **An MLST component for [bifrost](https://github.com/ssi-dk/bifrost) utilizing ariba**

Category | Value
--- | ---
name | ariba_mlst
version | 1.1
target | sample
type | pipeline
recommendation | recommended

### **Description**:
This preforms read based mapping against the mlst DB's. The mlst DB's are set via the species 
table in the mongoDB. This was previously a part of Analyzer, but is better maintained by 
splitting that component apart
#---------------------------------------------------------------------------------------------------

#-Options-------------------------------------------------------------------------------------------
# None
#---------------------------------------------------------------------------------------------------

#-Required resource files---------------------------------------------------------------------------
# Relative files are relative to this components folder if not using an absolute path, a recommended
# strategy for shrared resources is to symlink the resource from the root resources folder into
# the component resources folder. This will prevent
mlst_database_path: "resources/ariba_mlst_db" 
#---------------------------------------------------------------------------------------------------

#-Requirements to run component---------------------------------------------------------------------
requirements:
  sample:
    properties:
      species:
#---------------------------------------------------------------------------------------------------

