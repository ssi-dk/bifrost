"bifrostlib.datahandling"
*************************


Module Contents
===============


Classes
-------

+------------+--------------------------------------------------------------------------------------------+
| "BifrostO  | For schema datatypes                                                                       |
| bjectData  |                                                                                            |
| Type"      |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "ObjectID" | For schema datatypes                                                                       |
+------------+--------------------------------------------------------------------------------------------+
| "BifrostO  |                                                                                            |
| bjectRefe  |                                                                                            |
| rence"     |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "BifrostO  |                                                                                            |
| bject"     |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "Category" |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "Sample"   |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "Run"      |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "Componen  |                                                                                            |
| t"         |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "Host"     |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "SampleCo  |                                                                                            |
| mponent"   |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "RunCompo  |                                                                                            |
| nent"      |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+


Functions
---------

+------------+--------------------------------------------------------------------------------------------+
| "load_sch  | loads BIFROST_SCHEMA from bifrost.jsonc which is the basis for objects                     |
| ema"() →   |                                                                                            |
| Dict       |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "get_sche  | Get a object schema from the BIFROST_SCHEMA                                                |
| ma_object  |                                                                                            |
| "(object_  |                                                                                            |
| type: str, |                                                                                            |
| schema_ve  |                                                                                            |
| rsion:     |                                                                                            |
| str) →     |                                                                                            |
| Dict       |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "get_sche  | Get a datatype schema from the BIFROST_SCHEMA                                              |
| ma_dataty  |                                                                                            |
| pes"(prog  |                                                                                            |
| ram_type:  |                                                                                            |
| str,       |                                                                                            |
| datatype:  |                                                                                            |
| str) →     |                                                                                            |
| Dict       |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "get_sche  | Get the reference schema of from the BIFROST_SCHEMA                                        |
| ma_refere  |                                                                                            |
| nce"(refe  |                                                                                            |
| rence_typ  |                                                                                            |
| e: str, s  |                                                                                            |
| chema_ver  |                                                                                            |
| sion: str) |                                                                                            |
| → Dict     |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+

bifrostlib.datahandling.BIFROST_SCHEMA

bifrostlib.datahandling.load_schema() -> Dict

   loads BIFROST_SCHEMA from bifrost.jsonc which is the basis for
   objects

   Other Parameters:
      BIFROST_SCHEMA (dict): GLOBAL storing the BIFROST_SCHEMA

   Returns:
      Dict: json formatted schema

bifrostlib.datahandling.get_schema_object(object_type: str, schema_version: str) -> Dict

   Get a object schema from the BIFROST_SCHEMA

   Note:
      With how it’s organized references and datatypes need to be
      included with object

   Args:
      object_type (str): object type based on available objects in
      schema schema_version (str): the version of the object you want
      to work with

   Other Parameters:
      BIFROST_SCHEMA (dict): GLOBAL storing the BIFROST_SCHEMA

   Returns:
      Dict: The object schema with datatypes schema and references
      schema it may need

bifrostlib.datahandling.get_schema_datatypes(program_type: str, datatype: str) -> Dict

   Get a datatype schema from the BIFROST_SCHEMA

   Args:
      reference_type (str): reference type based on available
      references to objects in schema

   Other Parameters:
      BIFROST_SCHEMA (dict): GLOBAL storing the BIFROST_SCHEMA

   Returns:
      Dict: The datatype schema

bifrostlib.datahandling.get_schema_reference(reference_type: str, schema_version: str) -> Dict

   Get the reference schema of from the BIFROST_SCHEMA

   Args:
      reference_type (str): [description] schema_version (str):
      [description]

   Returns:
      Dict: [description]

class bifrostlib.datahandling.BifrostObjectDataType(program_type: str, datatype: str, _json: Dict)

   For schema datatypes

   __repr__(self)

      Return repr(self).

   __getitem__(self, key)

   __setitem__(self, key, value)

   __delitem__(self, key)

class bifrostlib.datahandling.ObjectID(_id: str)

   Bases: "bifrostlib.datahandling.BifrostObjectDataType"

   For schema datatypes

class bifrostlib.datahandling.BifrostObjectReference(reference_type: str, value: Dict, schema_version: str)

   reference_type

   __repr__(self)

      Return repr(self).

   __getitem__(self, key)

   __setitem__(self, key, value)

   __delitem__(self, key)

class bifrostlib.datahandling.BifrostObject(object_type: str, value: Dict, schema_version: str)

   json

   __repr__(self)

      Return repr(self).

   __getitem__(self, key)

   __setitem__(self, key, value)

   __delitem__(self, key)

   load(self, _id: ObjectID)

   load_from_reference(self, reference: BifrostObjectReference)

   save(self)

   delete(self)

   to_reference(self, additional_requirements: Dict = {})

class bifrostlib.datahandling.Category(schema_version='2.1')

   Bases: "bifrostlib.datahandling.BifrostObject"

class bifrostlib.datahandling.Sample(name: str = None, schema_version='2.1')

   Bases: "bifrostlib.datahandling.BifrostObject"

   properties :Category

   components :List[BifrostObjectReference]

   get_category(self, key: str)

   set_category(self, category: Category)

class bifrostlib.datahandling.Run(name: str = None, schema_version='2.1')

   Bases: "bifrostlib.datahandling.BifrostObject"

   samples :List[BifrostObjectReference]

   components :List[BifrostObjectReference]

   hosts :List[BifrostObjectReference]

class bifrostlib.datahandling.Component(name: str = None, schema_version='2.1')

   Bases: "bifrostlib.datahandling.BifrostObject"

class bifrostlib.datahandling.Host(name: str = None, schema_version='2.1')

   Bases: "bifrostlib.datahandling.BifrostObject"

class bifrostlib.datahandling.SampleComponent(sample_ref: BifrostObjectReference, component_ref: BifrostObjectReference, schema_version='2.1')

   Bases: "bifrostlib.datahandling.BifrostObject"

class bifrostlib.datahandling.RunComponent(run_ref: BifrostObjectReference, component_ref: BifrostObjectReference, schema_version='2.1')

   Bases: "bifrostlib.datahandling.BifrostObject"
