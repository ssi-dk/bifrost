import mongomock
import pymongo
from bifrostlib import datahandling
import datetime
import time


def test_bifrostlib():
    assert datahandling.test() == u"This is the bifrost lib"


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_post_sample():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    sample["_id"] = sample_db["_id"]
    assert sample == sample_db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_get_sample():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    sample["_id"] = sample_db["_id"]
    sample_received = datahandling.get_samples([str(sample["_id"])])[0]
    assert sample_received == sample


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_get_sample_using_component_id():
    component = {"name": "assemblatron"}
    component_db = datahandling.post_component(component)

    sample = {
        "name": "test_sample",
        "components": [
            component_db
        ]
    }
    sample_db = datahandling.post_sample(sample)
    sample["_id"] = sample_db["_id"]
    sample_received = datahandling.get_samples(
        component_ids=[str(component_db["_id"])])[0]
    assert sample_received == sample



@mongomock.patch(('mongodb://server.example.com:27017'))
def test_post_component():
    component = {"name": "assemblatron"}
    component_db = datahandling.post_component(component)
    component["_id"] = component_db["_id"]
    assert component == component_db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_delete_component():
    component = {"name": "assemblatron"}
    component_db = datahandling.post_component(component)
    
    deleted = datahandling.delete_component(str(component_db["_id"]))
    assert deleted == 1
# There is not "get component" without db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_post_sample_component():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    component = {"name": "assemblatron"}
    component_db = datahandling.post_component(component)
    s_c = {
        "sample": {"_id": sample_db["_id"]},
        "component": {"_id": component_db["_id"]},
        "results": {},
        "summary": {}
    }
    s_c_db = datahandling.post_sample_component(s_c)
    s_c["_id"] = s_c_db["_id"]
    assert s_c == s_c_db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_get_sample_components():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    component = {"name": "assemblatron"}
    component_db = datahandling.post_component(component)
    s_c = {
        "sample": {"_id": sample_db["_id"]},
        "component": {
            "_id": component_db["_id"],
            "name": component_db["name"]
        },
        "results": {},
        "summary": {},
        "setup_date": datetime.datetime.now()
    }
    s_c_db = datahandling.post_sample_component(s_c)
    # Testing this.
    s_c_2 = datahandling.get_sample_components(
        sample_ids=[str(sample_db["_id"])],
        component_names=["assemblatron"])
    assert s_c_2[0] == s_c_db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_delete_sample():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    component = {"name": "assemblatron"}
    component_db = datahandling.post_component(component)
    s_c = {
        "sample": {"_id": sample_db["_id"]},
        "component": {
            "_id": component_db["_id"],
            "name": component_db["name"]
        },
        "results": {},
        "summary": {},
        "setup_date": datetime.datetime.now()
    }
    datahandling.post_sample_component(s_c)

    # Testing this
    deleted = datahandling.delete_sample(str(sample_db["_id"]))
    assert deleted == 1


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_post_run():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    component = {"name": "assemblatron"}
    component_db = datahandling.post_component(component)
    s_c = {
        "sample": {"_id": sample_db["_id"]},
        "component": {
            "_id": component_db["_id"],
            "name": component_db["name"]
        },
        "results": {},
        "summary": {},
        "setup_date": datetime.datetime.now()
    }
    datahandling.post_sample_component(s_c)

    run = {
        "name": "test_run",
        "samples": [
            {"_id": sample_db["_id"]}
        ],
        "components": [
            {"_id": component_db["_id"]}
            # Missing name
        ]
    }
    run_db = datahandling.post_run(run)

    run["_id"] = run_db["_id"]
    
    assert run == run_db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_get_run():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    component = {"name": "assemblatron"}
    component_db = datahandling.post_component(component)
    s_c = {
        "sample": {"_id": sample_db["_id"]},
        "component": {
            "_id": component_db["_id"],
            "name": component_db["name"]
        },
        "results": {},
        "summary": {},
        "setup_date": datetime.datetime.now()
    }
    datahandling.post_sample_component(s_c)

    run = {
        "name": "test_run",
        "samples": [
            {"_id": sample_db["_id"]}
        ],
        "components": [
            {"_id": component_db["_id"]}
            # Missing name
        ]
    }
    run_db = datahandling.post_run(run)

    run_get = datahandling.get_runs(run_id=str(run_db["_id"]))[0]

    assert run_get == run_db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_delete_run():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    component = {"name": "assemblatron"}
    component_db = datahandling.post_component(component)
    s_c = {
        "sample": {"_id": sample_db["_id"]},
        "component": {
            "_id": component_db["_id"],
            "name": component_db["name"]
        },
        "results": {},
        "summary": {},
        "setup_date": datetime.datetime.now()
    }
    datahandling.post_sample_component(s_c)

    run = {
        "name": "test_run",
        "samples": [
            {"_id": sample_db["_id"]}
        ],
        "components": [
            {"_id": component_db["_id"]}
            # Missing name
        ]
    }
    run_db = datahandling.post_run(run)

    deleted = datahandling.delete_run(run_id=str(run_db["_id"]))

    assert deleted == 1


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_export_run():
    component = {"name": "assemblatron"}
    component_db = datahandling.post_component(component)
    component_2 = {"name": "whats_my_species"}
    component_db_2 = datahandling.post_component(component_2)
    sample = {
        "name": "test_sample",
        "components": [
            {"_id": component_db["_id"]},
            {"_id": component_db_2["_id"]}
        ]
    }
    sample_db = datahandling.post_sample(sample)
    s_c = {
        "sample": {"_id": sample_db["_id"]},
        "component": {
            "_id": component_db["_id"],
            "name": component_db["name"]
        },
        "results": {},
        "summary": {},
        "setup_date": datetime.datetime.now()
    }
    s_c_db = datahandling.post_sample_component(s_c)
    time.sleep(0.2)
    # we need to wait otherwise they share the same setup_date and
    # it screws up the order
    s_c_2 = {
        "sample": {"_id": sample_db["_id"]},
        "component": {
            "_id": component_db_2["_id"],
            "name": component_db_2["name"]
        },
        "results": {},
        "summary": {},
        "setup_date": datetime.datetime.now()
    }
    s_c_db_2 = datahandling.post_sample_component(s_c_2)

    run = {
        "name": "test_run",
        "samples": [
            {"_id": sample_db["_id"]}
        ],
        "components": [
            {"_id": component_db["_id"]},
            {"_id": component_db_2["_id"]}
            # Missing name
        ]
    }
    run_db = datahandling.post_run(run)

    run_export = datahandling.get_run_export(run_ids=[str(run_db["_id"])])

    run_id = str(run_db["_id"])

    run_expected = {
        run_id: {
            "components": [
                component_db_2,
                component_db
                ],
            "samples": [sample_db],
            "sample_components": [
                s_c_db_2,
                s_c_db
            ],
            "runs": [run_db]
        }
    }

    assert run_export == run_expected





@mongomock.patch(('mongodb://server.example.com:27017'))
def test_inport_run():
    component = {"name": "assemblatron"}
    component_db = datahandling.post_component(component)
    component_2 = {"name": "whats_my_species"}
    component_db_2 = datahandling.post_component(component_2)
    sample = {
        "name": "test_sample",
        "components": [
            {"_id": component_db["_id"]},
            {"_id": component_db_2["_id"]}
        ]
    }
    sample_db = datahandling.post_sample(sample)
    s_c = {
        "sample": {"_id": sample_db["_id"]},
        "component": {
            "_id": component_db["_id"],
            "name": component_db["name"]
        },
        "results": {},
        "summary": {},
        "setup_date": datetime.datetime.now()
    }
    s_c_db = datahandling.post_sample_component(s_c)

    time.sleep(0.2)
    # we need to wait otherwise they share the same setup_date and
    # it screws up the order
    
    s_c_2 = {
        "sample": {"_id": sample_db["_id"]},
        "component": {
            "_id": component_db_2["_id"],
            "name": component_db_2["name"]
        },
        "results": {},
        "summary": {},
        "setup_date": datetime.datetime.now()
    }
    s_c_db_2 = datahandling.post_sample_component(s_c_2)

    run = {
        "name": "test_run",
        "samples": [
            {"_id": sample_db["_id"]}
        ],
        "components": [
            {"_id": component_db["_id"]},
            {"_id": component_db_2["_id"]}
            # Missing name
        ]
    }
    run_db = datahandling.post_run(run)

    run_export = datahandling.get_run_export(run_ids=[str(run_db["_id"])])

    #drop database
    conn = pymongo.MongoClient('mongodb://server.example.com:27017')
    conn.drop_database("serumqc_prod")
    conn.close()

    imported = datahandling.post_run_export(run_export)

    get_items = {
        "samples": [],
        "sample_components": [],
        "runs": [],
        "components": []
    }

    run_dict = run_export[str(run_db["_id"])]

    for s in run_dict["samples"]:
        get_s = datahandling.get_samples(sample_ids=[str(s["_id"])])[0]
        get_items["samples"].append(get_s)

    for s_c in run_dict["sample_components"]:
        get_s_c = datahandling.get_sample_components(
            sample_component_ids=[str(s_c["_id"])])[0]
        get_items["sample_components"].append(get_s_c)

    for r in run_dict["runs"]:
        get_r = datahandling.get_runs(run_id=str(r["_id"]))[0]
        get_items["runs"].append(get_r)

    for c in run_dict["components"]:
        get_c = datahandling.get_components(component_ids=[str(c["_id"])])[0]
        get_items["components"].append(get_c)
    
    assert get_items == run_dict
