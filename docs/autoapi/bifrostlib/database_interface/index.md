# `bifrostlib.database_interface`

## Module Contents

### Functions

| `date_now`() → datetime.datetime

 | Get the current time as a datetime

 |
| `get_connection`() → pymongo.mongo_client.MongoClient

 | Get a connection to the DB

         |
| `close_connection`()

                                  | Closes DB connection

               |
| `pluralize`(name: str) → str

                          | Turns a string into a pluralized form. For example sample -> samples and property -> properties

 |
| `json_to_bson`(json_object: Dict) → Dict

              | Converts a json dict to bson dict

                                                               |
| `bson_to_json`(bson_object: Dict) → Dict

              | Converts a bson dict to json dict

                                                               |
| `load`(object_type: str, _id: Dict) → Dict

            | Loads an object based on it’s id from the DB

                                                    |
| `save`(object_type, object_value: Dict) → Dict

        | Saves a object to the DB

                                                                        |
| `delete`(object_type, _id: Dict) → bool

               | Deletes a object from the DB based on it’s id

                                                   |

### bifrostlib.database_interface.date_now()
Get the current time as a datetime

**NOTE**: Needed to keep the same date in python and mongo, as mongo rounds to millisecond


* **Returns**

    Current time rounded to miliseconds



### bifrostlib.database_interface.CONNECTION()

### bifrostlib.database_interface.get_connection()
Get a connection to the DB


* **Other Parameters**

    
    * **CONNECTION** (*pymongo.MongoClient*) – GLOBAL storing connection


    * **BIFROST_DB_KEY** – (ENV) This is a environmental variable taken from the system



* **Returns**

    Sets global CONNECTION string based on env var BIFROST_DB_KEY



* **Return type**

    pymongo.mongo_client.MongoClient



* **Raises**

    [**ValueError**](https://docs.python.org/3/library/exceptions.html#ValueError) – If DB is not set properly



### bifrostlib.database_interface.close_connection()
Closes DB connection


* **Other Parameters**

    **CONNECTION** (*pymongo.MongoClient*) – GLOBAL storing connection



### bifrostlib.database_interface.pluralize(name: [str](https://docs.python.org/3/library/stdtypes.html#str))
Turns a string into a pluralized form. For example sample -> samples and property -> properties


* **Parameters**

    **name** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – A non plural string to turn into it’s plural



* **Returns**

    The pluralized form of the string.



* **Return type**

    [str](https://docs.python.org/3/library/stdtypes.html#str)



### bifrostlib.database_interface.json_to_bson(json_object: Dict)
Converts a json dict to bson dict


* **Parameters**

    **json_object** (*Dict*) – A json formatted dict



* **Returns**

    A bson formatted dict



* **Return type**

    Dict



### bifrostlib.database_interface.bson_to_json(bson_object: Dict)
Converts a bson dict to json dict


* **Parameters**

    **bson_object** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – A bson formatted dict



* **Returns**

    A json formatted dict



* **Return type**

    Dict



### bifrostlib.database_interface.load(object_type: [str](https://docs.python.org/3/library/stdtypes.html#str), _id: Dict)
Loads an object based on it’s id from the DB

**NOTE**: Inputs and outputs are json dict but database works on bson dicts


* **Parameters**

    
    * **object_type** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – A bifrost object type found in the database as a collection


    * **_id** (*Dict*) – json formatted objectid {“$oid”: <value>}



* **Returns**

    json formatted dict of the object



* **Return type**

    Dict



* **Raises**

    [**AssertionError**](https://docs.python.org/3/library/exceptions.html#AssertionError) – If db contains a duplicate _id



### bifrostlib.database_interface.save(object_type, object_value: Dict)
Saves a object to the DB

**NOTE**: Inputs and outputs are json dict but database works on bson dicts


* **Parameters**

    
    * **object_type** – A bifrost object type found in the database as a collection


    * **object_value** – json formatted object



* **Returns**

    json formatted dict of the object with objectid



* **Return type**

    Dict



* **Raises**

    [**KeyError**](https://docs.python.org/3/library/exceptions.html#KeyError) – If object_type not in DB



### bifrostlib.database_interface.delete(object_type, _id: Dict)
Deletes a object from the DB based on it’s id

**NOTE**: This only removes the objects and doesn’t handle dependencies or dangling objects


* **Parameters**

    
    * **object_type** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – A bifrost object type found in the database as a collection


    * **_id** (*Dict*) – json formatted objectid {“$oid”: <value>}



* **Returns**

    Successfully deleted | Failure to delete



* **Return type**

    [bool](https://docs.python.org/3/library/functions.html#bool)



* **Raises**

    [**KeyError**](https://docs.python.org/3/library/exceptions.html#KeyError) – If object_type not in DB
