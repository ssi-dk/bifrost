# `bifrostlib.datahandling`

## Module Contents

### Classes

| `BifrostObjectDataType`
 | For schema datatypes

 |
| `ObjectID`
                                            | For schema datatypes

                                                                            |
| `BifrostObjectReference`
                              | 

                                                                                                |
| `BifrostObject`
                                       | 

                                                                                                |
| `Category`
                                            | 

                                                                                                |
| `Sample`
                                              | 

                                                                                                |
| `Run`
                                                 | 

                                                                                                |
| `Component`
                                           | 

                                                                                                |
| `Host`
                                                | 

                                                                                                |
| `SampleComponent`
                                     | 

                                                                                                |
| `RunComponent`
                                        | 

                                                                                                |
### Functions

| `load_schema`() → Dict

                                | loads BIFROST_SCHEMA from bifrost.jsonc which is the basis for objects

                          |
| `get_schema_object`(object_type: str, schema_version: str) → Dict

 | Get a object schema from the BIFROST_SCHEMA

                                                     |
| `get_schema_datatypes`(program_type: str, datatype: str) → Dict

   | Get a datatype schema from the BIFROST_SCHEMA

                                                   |
| `get_schema_reference`(reference_type: str, schema_version: str) → Dict

 | Get the reference schema of from the BIFROST_SCHEMA

                                             |

### bifrostlib.datahandling.BIFROST_SCHEMA()

### bifrostlib.datahandling.load_schema()
loads BIFROST_SCHEMA from bifrost.jsonc which is the basis for objects


* **Other Parameters**

    **BIFROST_SCHEMA** (*dict*) – GLOBAL storing the BIFROST_SCHEMA



* **Returns**

    json formatted schema



* **Return type**

    Dict



### bifrostlib.datahandling.get_schema_object(object_type: [str](https://docs.python.org/3/library/stdtypes.html#str), schema_version: [str](https://docs.python.org/3/library/stdtypes.html#str))
Get a object schema from the BIFROST_SCHEMA

**NOTE**: With how it’s organized references and datatypes need to be included with object


* **Parameters**

    
    * **object_type** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – object type based on available objects in schema


    * **schema_version** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – the version of the object you want to work with



* **Other Parameters**

    **BIFROST_SCHEMA** (*dict*) – GLOBAL storing the BIFROST_SCHEMA



* **Returns**

    The object schema with datatypes schema and references schema it may need



* **Return type**

    Dict



### bifrostlib.datahandling.get_schema_datatypes(program_type: [str](https://docs.python.org/3/library/stdtypes.html#str), datatype: [str](https://docs.python.org/3/library/stdtypes.html#str))
Get a datatype schema from the BIFROST_SCHEMA


* **Parameters**

    **reference_type** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – reference type based on available references to objects in schema



* **Other Parameters**

    **BIFROST_SCHEMA** (*dict*) – GLOBAL storing the BIFROST_SCHEMA



* **Returns**

    The datatype schema



* **Return type**

    Dict



### bifrostlib.datahandling.get_schema_reference(reference_type: [str](https://docs.python.org/3/library/stdtypes.html#str), schema_version: [str](https://docs.python.org/3/library/stdtypes.html#str))
Get the reference schema of from the BIFROST_SCHEMA


* **Parameters**

    
    * **reference_type** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – [description]


    * **schema_version** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – [description]



* **Returns**

    [description]



* **Return type**

    Dict



### class bifrostlib.datahandling.BifrostObjectDataType(program_type: [str](https://docs.python.org/3/library/stdtypes.html#str), datatype: [str](https://docs.python.org/3/library/stdtypes.html#str), _json: Dict)
For schema datatypes


#### \__repr__(self)
Return repr(self).


#### \__getitem__(self, key)

#### \__setitem__(self, key, value)

#### \__delitem__(self, key)

### class bifrostlib.datahandling.ObjectID(_id: [str](https://docs.python.org/3/library/stdtypes.html#str))
Bases: `bifrostlib.datahandling.BifrostObjectDataType`

For schema datatypes


### class bifrostlib.datahandling.BifrostObjectReference(reference_type: [str](https://docs.python.org/3/library/stdtypes.html#str), value: Dict, schema_version: [str](https://docs.python.org/3/library/stdtypes.html#str))

#### \__repr__(self)
Return repr(self).


#### \__getitem__(self, key)

#### \__setitem__(self, key, value)

#### \__delitem__(self, key)

#### property reference_type(self)

### class bifrostlib.datahandling.BifrostObject(object_type: [str](https://docs.python.org/3/library/stdtypes.html#str), value: Dict, schema_version: [str](https://docs.python.org/3/library/stdtypes.html#str))

#### \__repr__(self)
Return repr(self).


#### \__getitem__(self, key)

#### \__setitem__(self, key, value)

#### \__delitem__(self, key)

#### property json(self)

#### load(self, _id: ObjectID)

#### load_from_reference(self, reference: BifrostObjectReference)

#### save(self)

#### delete(self)

#### to_reference(self, additional_requirements: Dict = {})

### class bifrostlib.datahandling.Category(schema_version='2.1')
Bases: `bifrostlib.datahandling.BifrostObject`


### class bifrostlib.datahandling.Sample(name: [str](https://docs.python.org/3/library/stdtypes.html#str) = None, schema_version='2.1')
Bases: `bifrostlib.datahandling.BifrostObject`


#### property properties(self)

#### property components(self)

#### get_category(self, key: [str](https://docs.python.org/3/library/stdtypes.html#str))

#### set_category(self, category: Category)

### class bifrostlib.datahandling.Run(name: [str](https://docs.python.org/3/library/stdtypes.html#str) = None, schema_version='2.1')
Bases: `bifrostlib.datahandling.BifrostObject`


#### property samples(self)

#### property components(self)

#### property hosts(self)

### class bifrostlib.datahandling.Component(name: [str](https://docs.python.org/3/library/stdtypes.html#str) = None, schema_version='2.1')
Bases: `bifrostlib.datahandling.BifrostObject`


### class bifrostlib.datahandling.Host(name: [str](https://docs.python.org/3/library/stdtypes.html#str) = None, schema_version='2.1')
Bases: `bifrostlib.datahandling.BifrostObject`


### class bifrostlib.datahandling.SampleComponent(sample_ref: BifrostObjectReference, component_ref: BifrostObjectReference, schema_version='2.1')
Bases: `bifrostlib.datahandling.BifrostObject`


### class bifrostlib.datahandling.RunComponent(run_ref: BifrostObjectReference, component_ref: BifrostObjectReference, schema_version='2.1')
Bases: `bifrostlib.datahandling.BifrostObject`
