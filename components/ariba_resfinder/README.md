#-Information---------------------------------------------------------------------------------------
name: ariba_resfinder
version: 1.1
target: sample
type: pipeline
recommendation: recommended
description: >
  This preforms read based mapping against the resfinder DB's.
#---------------------------------------------------------------------------------------------------

#-Options-------------------------------------------------------------------------------------------
# -None
#---------------------------------------------------------------------------------------------------

#-Required resource files---------------------------------------------------------------------------
# Relative files are relative to this components folder if not using an absolute path
abricate_resfinder_database: "resources/abricate_resfinder_db" 
ariba_resfinder_database: "resources/ariba_resfinder_db"
#---------------------------------------------------------------------------------------------------

#-Requirements to run component---------------------------------------------------------------------
requirements:
  # None
#---------------------------------------------------------------------------------------------------