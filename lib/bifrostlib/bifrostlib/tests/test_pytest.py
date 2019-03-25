import mongomock
from bifrostlib import datahandling
import datetime


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
    component_db = datahandling.save_component_to_db(component)

    sample = {
        "name": "test_sample",
        "components": [
            component_db
        ]
    }
    sample_db = datahandling.post_sample(sample)
    sample["_id"] = sample_db["_id"]
    sample_received = datahandling.get_samples(
        [str(component_db["_id"])])[0]
    assert sample_received == sample



@mongomock.patch(('mongodb://server.example.com:27017'))
def test_post_component():
    component = {"name": "assemblatron"}
    component_db = datahandling.save_component_to_db(component)
    component["_id"] = component_db["_id"]
    assert component == component_db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_delete_component():
    component = {"name": "assemblatron"}
    component_db = datahandling.save_component_to_db(component)
    
    deleted = datahandling.delete_component(str(component_db["_id"]))
    assert deleted == 1
# There is not "get component" without db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_post_sample_component():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    component = {"name": "assemblatron"}
    component_db = datahandling.save_component_to_db(component)
    s_c = {
        "sample": {"_id": sample_db["_id"]},
        "component": {"_id": component_db["_id"]},
        "results": {},
        "summary": {}
    }
    s_c_db = datahandling.save_sample_component_to_db(s_c)
    s_c["_id"] = s_c_db["_id"]
    assert s_c == s_c_db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_get_sample_components():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    component = {"name": "assemblatron"}
    component_db = datahandling.save_component_to_db(component)
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
    s_c_db = datahandling.save_sample_component_to_db(s_c)
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
    component_db = datahandling.save_component_to_db(component)
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
    datahandling.save_sample_component_to_db(s_c)

    # Testing this
    deleted = datahandling.delete_sample(str(sample_db["_id"]))
    assert deleted == 1


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_post_run():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    component = {"name": "assemblatron"}
    component_db = datahandling.save_component_to_db(component)
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
    datahandling.save_sample_component_to_db(s_c)

    run = {
        "name": "test_run",
        "samples": [
            "_id": sample_db["_id"]
        ]
        "components": [
            "_id": component_db["_id"]
            # Missing name
        ]
    }
    run_db = datahandling.post_run(run)

    run["_id"] == run_db["_id"]
    
    assert run == run_db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_get_run():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    component = {"name": "assemblatron"}
    component_db = datahandling.save_component_to_db(component)
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
    datahandling.save_sample_component_to_db(s_c)

    run = {
        "name": "test_run",
        "samples": [
            "_id": sample_db["_id"]
        ]
        "components": [
            "_id": component_db["_id"]
            # Missing name
        ]
    }
    run_db = datahandling.post_run(run)

    run_get = datahandling.get_runs(run_id=run_db["_id"])[0]

    assert run_get == run_db


@mongomock.patch(('mongodb://server.example.com:27017'))
def test_delete_run():
    sample = {"name": "test_sample"}
    sample_db = datahandling.post_sample(sample)
    component = {"name": "assemblatron"}
    component_db = datahandling.save_component_to_db(component)
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
    datahandling.save_sample_component_to_db(s_c)

    run = {
        "name": "test_run",
        "samples": [
            "_id": sample_db["_id"]
        ]
        "components": [
            "_id": component_db["_id"]
            # Missing name
        ]
    }
    run_db = datahandling.post_run(run)

    deleted = datahandling.delete_run(run_id=run_db["_id"])[0]

    assert deleted == 1
