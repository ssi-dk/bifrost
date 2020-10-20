"bifrostlib.database_interface"
*******************************


Module Contents
===============


Functions
---------

+------------+--------------------------------------------------------------------------------------------+
| "date_now  | Get the current time as a datetime                                                         |
| "() → dat  |                                                                                            |
| etime.dat  |                                                                                            |
| etime      |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "get_conn  | Get a connection to the DB                                                                 |
| ection"()  |                                                                                            |
| → pymongo  |                                                                                            |
| .mongo_cl  |                                                                                            |
| ient.Mong  |                                                                                            |
| oClient    |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "close_co  | Closes DB connection                                                                       |
| nnection"  |                                                                                            |
| ()         |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "pluraliz  | Turns a string into a pluralized form. For example sample -> samples and property ->       |
| e"(name:   | properties                                                                                 |
| str) → str |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "json_to_  | Converts a json dict to bson dict                                                          |
| bson"(jso  |                                                                                            |
| n_object:  |                                                                                            |
| Dict) →    |                                                                                            |
| Dict       |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "bson_to_  | Converts a bson dict to json dict                                                          |
| json"(bso  |                                                                                            |
| n_object:  |                                                                                            |
| Dict) →    |                                                                                            |
| Dict       |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "load"(ob  | Loads an object based on it’s id from the DB                                               |
| ject_type: |                                                                                            |
| str, _id:  |                                                                                            |
| Dict) →    |                                                                                            |
| Dict       |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "save"(ob  | Saves a object to the DB                                                                   |
| ject_type, |                                                                                            |
| object_va  |                                                                                            |
| lue: Dict) |                                                                                            |
| → Dict     |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+
| "delete"(  | Deletes a object from the DB based on it’s id                                              |
| object_ty  |                                                                                            |
| pe, _id:   |                                                                                            |
| Dict) →    |                                                                                            |
| bool       |                                                                                            |
+------------+--------------------------------------------------------------------------------------------+

bifrostlib.database_interface.date_now() -> datetime.datetime

   Get the current time as a datetime

   Note:
      Needed to keep the same date in python and mongo, as mongo
      rounds to millisecond

   Returns:
      Current time rounded to miliseconds

bifrostlib.database_interface.CONNECTION

bifrostlib.database_interface.get_connection() -> pymongo.mongo_client.MongoClient

   Get a connection to the DB

   Other Parameters:
      CONNECTION (pymongo.MongoClient): GLOBAL storing connection
      BIFROST_DB_KEY: (ENV) This is a environmental variable taken
      from the system

   Returns:
      pymongo.mongo_client.MongoClient: Sets global CONNECTION string
      based on env var BIFROST_DB_KEY

   Raises:
      ValueError: If DB is not set properly

bifrostlib.database_interface.close_connection()

   Closes DB connection

   Other Parameters:
      CONNECTION (pymongo.MongoClient): GLOBAL storing connection

bifrostlib.database_interface.pluralize(name: str) -> str

   Turns a string into a pluralized form. For example sample ->
   samples and property -> properties

   Args:
      name (str): A non plural string to turn into it’s plural

   Returns:
      str: The pluralized form of the string.

bifrostlib.database_interface.json_to_bson(json_object: Dict) -> Dict

   Converts a json dict to bson dict

   Args:
      json_object (Dict): A json formatted dict

   Returns:
      Dict: A bson formatted dict

bifrostlib.database_interface.bson_to_json(bson_object: Dict) -> Dict

   Converts a bson dict to json dict

   Args:
      bson_object (str): A bson formatted dict

   Returns:
      Dict: A json formatted dict

bifrostlib.database_interface.load(object_type: str, _id: Dict) -> Dict

   Loads an object based on it’s id from the DB

   Note:
      Inputs and outputs are json dict but database works on bson
      dicts

   Args:
      object_type (str): A bifrost object type found in the database
      as a collection _id (Dict): json formatted objectid {“$oid”:
      <value>}

   Returns:
      Dict: json formatted dict of the object

   Raises:
      AssertionError: If db contains a duplicate _id

bifrostlib.database_interface.save(object_type, object_value: Dict) -> Dict

   Saves a object to the DB

   Note:
      Inputs and outputs are json dict but database works on bson
      dicts

   Args:
      object_type: A bifrost object type found in the database as a
      collection object_value: json formatted object

   Returns:
      Dict: json formatted dict of the object with objectid

   Raises:
      KeyError: If object_type not in DB

bifrostlib.database_interface.delete(object_type, _id: Dict) -> bool

   Deletes a object from the DB based on it’s id

   Note:
      This only removes the objects and doesn’t handle dependencies or
      dangling objects

   Args:
      object_type (str): A bifrost object type found in the database
      as a collection _id (Dict): json formatted objectid {“$oid”:
      <value>}

   Returns:
      bool: Successfully deleted | Failure to delete

   Raises:
      KeyError: If object_type not in DB
