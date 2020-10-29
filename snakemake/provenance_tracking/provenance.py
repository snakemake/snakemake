__author__ = "Alban Gaignard"
__email__ = "alban.gaignard@univ-nantes.fr"

import uuid
import getpass
import os
import hashlib
import prov.model as prov
import datetime

from snakemake.logging import logger


class PROVMgr(object):
    __instance = None
    __is_active = False

    whoami = getpass.getuser()

    def __new__(self):
        if PROVMgr.__instance is None:
            PROVMgr.__instance = object.__new__(self)
        return PROVMgr.__instance

    def __init__(self):
        self.document = prov.ProvDocument()
        self.document.set_default_namespace("http://snakemake-provenance#")
        self.document.add_namespace("xsd", "http://www.w3.org/2001/XMLSchema#")
        self.document.add_namespace(
            "crypto",
            "http://id.loc.gov/vocabulary/preservation/cryptographicHashFunctions#",
        )
        self.document.add_namespace("rdfs", "http://www.w3.org/2000/01/rdf-schema#")
        self.document.add_namespace("prov", "http://www.w3.org/ns/prov#")

        ag = self.document.agent(self.whoami)

        self.bundle = self.document.bundle(self.gen_URI())
        self.bundle.wasAttributedTo(agent=ag.identifier, entity=self.bundle.identifier)
        self.wfexec = self.bundle.activity(self.gen_URI(), startTime=self.__get_iso_now)

    def set_activate(self, is_active):
        self.__is_active = is_active

    def is_active(self):
        return self.__is_active

    def add_activity(self, tool_name, job_uri, input_id_list, cmd):
        activity = self.bundle.activity(
            identifier=str(job_uri), startTime=self.__get_iso_now
        )
        activity.wasAssociatedWith(agent=tool_name)

        if cmd:
            activity.add_attributes({"rdfs:comment": cmd})

        for in_id in input_id_list:
            input = self.bundle.entity(in_id)
            try:
                if os.path.isfile(os.path.abspath(in_id)):
                    digest = hashlib.sha512(
                        open(os.path.abspath(in_id), "rb").read()
                    ).hexdigest()
                    input.add_attributes({"crypto:sha512": digest})
            except Exception as e:
                logger.error(
                    "Error while computing the SHA512 fingerprint for " + str(in_id)
                )
            activity.used(entity=input)

    def add_output(self, output_id, input_id_list, tool_name, job_uri):
        output = self.bundle.entity(output_id)
        output.wasGeneratedBy(activity=str(job_uri))
        output.wasAttributedTo(agent=str(tool_name))
        # rdfs:label ?

        try:
            if os.path.isfile(os.path.abspath(output_id)):
                digest = hashlib.sha512(
                    open(os.path.abspath(output_id), "rb").read()
                ).hexdigest()
                output.add_attributes({"crypto:sha512": digest})
        except Exception as e:
            logger.error(
                "Error while computing the SHA512 fingerprint for " + str(output_id)
            )

        for in_id in input_id_list:
            output.wasDerivedFrom(usedEntity=in_id)

    def terminate_wf_exec(self):
        self.wfexec.set_time(
            startTime=self.wfexec.get_startTime(), endTime=datetime.datetime.now()
        )
        self.document.serialize(destination="provenance.trig", format="rdf")
        self.document.serialize(destination="provenance.json", format="json")

    @property
    def __get_iso_now(self):
        return datetime.datetime.now().isoformat()

    @property
    def __get_xsd_now(self):
        return datetime.datetime.now().isoformat() + '"^^xsd:dateTime'

    def gen_URI(self, prefix="", name=""):
        if str(prefix + name) == "":
            return str(uuid.uuid4())
        else:
            return prefix + name + "-" + str(uuid.uuid4())


############ Singleton class instanciation ###########
provenance_manager = PROVMgr()
