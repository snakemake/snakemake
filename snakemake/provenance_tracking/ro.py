__author__ = "Alban Gaignard, Johannes Köster"
__copyright__ = "Copyright 2020, Alban Gaignard, Johannes Köster"
__license__ = "MIT"

# from snakemake.logging import logger
import json


class RO:

    root = {
        "@type": "CreativeWork",
        "@id": "ro-crate-metadata.jsonld",
        "identifier": "ro-crate-metadata.jsonld",
        "about": {"@id": "./"},
    }

    def __init__(self, name):
        self.name = name

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)


if __name__ == "__main__":
    ro1 = RO(name="toto")
    print(ro1.toJSON())
    print(ro1.root)
