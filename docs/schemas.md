# Data Structure

- [Data Structure](#data-structure)
  - [Overview](#overview)
  - [Object references](#object-references)
  - [Objects](#objects)
    - [**sample** (genomicsamples)](#sample-genomicsamples)
    - [**run** (collections)](#run-collections)
    - [**component** (pipelines)](#component-pipelines)
    - [**properties/category** (categories)](#propertiescategory-categories)
    - [**sample_components** (results_of_pipeline_on_genomicsample)](#sample_components-results_of_pipeline_on_genomicsample)
    - [**epis** (hosts)](#epis-hosts)
    - [** epievents** (events)](#-epievents-events)
    - [**run_components** (results_of_pipeline_on_collection)](#run_components-results_of_pipeline_on_collection)

## Overview
The bifrost DB is key to organization of genomic/epi data and keeping track of how all of it is processed. Currently things are built around the **sample** object. You'll see the headings below include a name in brackets, this is a more descriptive name for each collection and future refactoring may reflect this.

## Object references
The following is a summary of how the objects interact with one another.
```
sample -> {[components], [categories], epi}
run -> {[samples], [components], [epis]}
component -> {[categories]}
category -> {}
sample_component -> {[categories], sample, component}
epi -> {[samples]}
run_component -> {[categories], run, component}
```
## Objects
These are the logical objects the DB is built around and often the Collection name in MongoDB

### **sample** (genomicsamples)
Stores the genomic data and is a base unit to save paired_read data and summarized results from the **components** run on them
```
[componentsObj]
[categoryObj]
{epiObj} # TODO
_dict = {
    _id: "ObjectID"
    name: "string"
    version: 
        schema: ["string"]
    metadata:
        created_at: "DateTime"
        updated_at: "DateTime"
    components: ["ComponentRefObj", "componentObj.status"]
    properties: {"CategoryRefObj", "categoryObj.summary"}
    report: {"CategoryRefObj", "categoryObj.report"}
    epi: {"EpiRefObj"} # TODO
}
```
NOTE: Epi hasn't been implented yet

### **run** (collections)
Organization for **samples** and **epi** including what **components** should run on the collection
```
[sampleObj]
[componentObj]
[epiObj]
_dict = {
    _id: "ObjectID"
    name: "string"
    version:
        schema: ["string"]
    metadata:
        created_at: "DateTime"
        updated_at: "DateTime"
    samples: ["SampleRefObj"]
    components: ["ComponentRefObj"]
    epis: ["EpiRefObj"] # TODO
    type: "string"
    path: "string
    comments: "string"
    issues:
        samples:
            duplicate_samples: ["string"]
            modified_samples: ["string"]
            unused_files: ["string"]
            samples_without_reads: ["string"]
            samples_without_metadata: ["string"]
}
```

### **component** (pipelines)
The storage of aditional information against a **run** or **sample**, this is typically a pipeline run against the data but can also include answers for questions against the DB (i.e. did the sample pass QC based off of other components run).
```
[categoryObj]
_dict = {
    _id: "ObjectID"
    name: "string"
    display_name: "string"
    version:
        schema: ["string"]
        code: "string"
        resources: "string"
    metadata:
        created_at: "DateTime"
        updated_at: "DateTime"
    category: ["CategoryRefObj"]
    db_values_changes:
        files:
        sample?:
            properties: {"CategoryRefObj", "categoryObj.summary"}
            report: {"CategoryRefObj", "categoryObj.report"}
        sample_component?:
            summary:  {"CategoryRefObj", "categoryObj.summary"}
            results: {}
        run?:
            properties: {"CategoryRefObj", "categoryObj.summary"}
            report: {"CategoryRefObj", "categoryObj.report"}
        run_component?:
            summary: {"CategoryRefObj", "categoryObj.summary"}
            results: {}
    details:
        target: "string" in ["sample", "run"]
        description: "string"
    options: {}
    resources: {}
    requirements:
        sample?: {}
        run?: {}
        component?: {}
    install:
        path: "string"
        dockerfile: "string"
}
```

### **properties/category** (categories)
A broader organization for **components** that do the same analysis but in different ways for example running two different assemblers against the same sample would have different results but cover the same category of denovo_assembly.
```
[categoryObj]
_dict = {
    name: "string" in ["paired_reads", "species_detection", "denovo_assembly", "stamper", "mlst", "plasmid", "resistance", "virulence", "size_check"]
    summary: {}
    report: {}
}
```
NOTE: This is a conceptual object and not a DB object

### **sample_components** (results_of_pipeline_on_genomicsample)
When running a **component** against a **sample** the results go here, the full results can be accessed but the summarized value and report also go into the sample **object**
```
{sampleObj}
{componentObj}
[categoryObj]
_dict = {
    _id: "ObjectID"
    version:
        schema: ["string"]
    metadata:
        created_at: "DateTime"
        updated_at: "DateTime"
    sample: {"SampleRefObj"}
    component: {"ComponentRefObj"}
    path: "string"
    status: "string" in ["Success", "Failure", "Initialized", "Requirements not met", "Queued", "Running"]
    properties: {"CategoryRef", "categoryObj.summary"}
    report: {"CategoryRef", "categoryObj.report"}
    results: {}
}
```

### **epis** (hosts)
This will be the object to store epidemiology data that references genomic **samples**
```
[sampleObj]
[eventObj]
_dict = {
    _id: "ObjectID"
    name: "string"
    version:
        schema: ["string"]
    metadata:
        created_at: "DateTime"
        updated_at: "DateTime"
    samples: ["SampleObjRef"]
    events: ["eventObj"]
    TBD # Need to work out more
}
```
NOTE: Not yet implemented, probably will need an eventObj 

### ** epievents** (events)
This will be the object to store the events related to an **epi**. It will likely not be a DB component.
```
_dict = {
    TBD # Need to work out more
}
```
NOTE: Not yet implemented

### **run_components** (results_of_pipeline_on_collection)
This will be the object which stores results of programs against a **run** (collection) of samples. An example of this would be SNP calling
```
{runObj}
{componentObj}
[categoryObj]
_dict = {
    _id: "ObjectID"
    version:
        schema: ["string"]
    metadata:
        created_at: "DateTime"
        updated_at: "DateTime"
    run: {"runRefObj"}
    component: {"ComponentRefObj"}
    path: "string"
    status: "string" in ["Success", "Failure", "Initialized", "Requirements not met", "Queued", "Running"]
    properties: {"CategoryRef", "categoryObj.summary"}
    report: {"CategoryRef", "categoryObj.report"}
    results: {}
}
```
NOTE: Not yet implemented