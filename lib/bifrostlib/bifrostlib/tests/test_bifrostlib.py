import mongomock
from bifrostlib import datahandling


def test_bifrostlib():
    assert datahandling.test() == u"This is the bifrost lib"


@mongomock.patch(('server.example.com', 27017))
def test_post_sample():
    sample = {"name": "test_sample"}
    sample_db = datahandling.save_sample_to_db(sample)
    sample["_id"] = sample_db["_id"]
    assert sample == sample_db
