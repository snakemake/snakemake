__author__ = "Christopher Tomkins-Tinch"
__copyright__ = "Copyright 2017, Christopher Tomkins-Tinch"
__email__ = "tomkinsc@broadinstitute.org"
__license__ = "MIT"

# built-ins
import time
import os
import re
import json
import logging
import xml.etree.ElementTree as ET

# module-specific
from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider
from snakemake.exceptions import WorkflowError, NCBIFileException
from snakemake.logging import logger

try:
    # third-party modules
    from Bio import Entrez
except ImportError as e:
    raise WorkflowError(
        "The Python package 'biopython' needs to be installed to use NCBI Entrez remote() file functionality. %s"
        % e.msg
    )


class RemoteProvider(AbstractRemoteProvider):
    def __init__(
        self,
        *args,
        keep_local=False,
        stay_on_remote=False,
        is_default=False,
        email=None,
        **kwargs
    ):
        super(RemoteProvider, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            is_default=is_default,
            email=email,
            **kwargs
        )
        self._ncbi = NCBIHelper(*args, email=email, **kwargs)

    def remote_interface(self):
        return self._ncbi

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "ncbi://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        return ["ncbi://"]

    def search(
        self, query, *args, db="nuccore", idtype="acc", retmode="json", **kwargs
    ):
        return list(
            self._ncbi.search(
                query, *args, db=db, idtype=idtype, retmode=retmode, **kwargs
            )
        )


class RemoteObject(AbstractRemoteObject):
    """ This is a class to interact with NCBI / GenBank.
    """

    def __init__(
        self,
        *args,
        keep_local=False,
        stay_on_remote=False,
        provider=None,
        email=None,
        db=None,
        rettype=None,
        retmode=None,
        **kwargs
    ):
        super(RemoteObject, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            provider=provider,
            email=email,
            db=db,
            rettype=rettype,
            retmode=retmode,
            **kwargs
        )
        if provider:
            self._ncbi = provider.remote_interface()
        else:
            self._ncbi = NCBIHelper(*args, email=email, **kwargs)

        if db and not self._ncbi.is_valid_db(db):
            raise NCBIFileException(
                "DB specified is not valid. Options include: {dbs}".format(
                    dbs=", ".join(self._ncbi.valid_dbs)
                )
            )
        else:
            self.db = db

        self.rettype = rettype
        self.retmode = retmode
        self.kwargs = kwargs

    # === Implementations of abstract class members ===

    def exists(self):
        if not self.retmode or not self.rettype:
            likely_request_options = self._ncbi.guess_db_options_for_extension(
                self.file_ext, db=self.db, rettype=self.rettype, retmode=self.retmode
            )
            self.db = likely_request_options["db"]
            self.retmode = likely_request_options["retmode"]
            self.rettype = likely_request_options["rettype"]
        return self._ncbi.exists(self.accession, db=self.db)

    def mtime(self):
        if self.exists():
            return self._ncbi.mtime(self.accession, db=self.db)
        else:
            raise NCBIFileException(
                "The record does not seem to exist remotely: %s" % self.accession
            )

    def size(self):
        if self.exists():
            return self._ncbi.size(self.accession, db=self.db)
        else:
            return self._iofile.size_local

    def download(self):
        if self.exists():
            self._ncbi.fetch_from_ncbi(
                [self.accession],
                os.path.dirname(self.accession),
                rettype=self.rettype,
                retmode=self.retmode,
                file_ext=self.file_ext,
                db=self.db,
                **self.kwargs
            )
        else:
            raise NCBIFileException(
                "The record does not seem to exist remotely: %s" % self.accession
            )

    def upload(self):
        raise NCBIFileException(
            "Upload is not permitted for the NCBI remote provider. Is an output set to NCBI.RemoteProvider.remote()?"
        )

    @property
    def list(self):
        raise NCBIFileException(
            "The NCBI Remote Provider does not currently support list-based operations like glob_wildcards()."
        )

    @property
    def accession(self):
        accession, version, file_ext = self._ncbi.parse_accession_str(self.local_file())
        return accession + "." + version

    @property
    def file_ext(self):
        accession, version, file_ext = self._ncbi.parse_accession_str(self.local_file())
        return file_ext

    @property
    def version(self):
        accession, version, file_ext = self._ncbi.parse_accession_str(self.local_file())
        return version


class NCBIHelper(object):
    def __init__(self, *args, email=None, **kwargs):
        if not email:
            raise NCBIFileException(
                "An e-mail address must be provided to either the remote file or the RemoteProvider() as email=<your_address>. The NCBI requires e-mail addresses for queries."
            )

        self.email = email
        self.entrez = Entrez
        self.entrez.email = self.email
        self.entrez.tool = "Snakemake"

        # valid NCBI Entrez efetch options
        # via https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
        self.efetch_options = {
            "bioproject": [{"rettype": "xml", "retmode": "xml", "ext": "xml"}],
            "biosample": [
                {"rettype": "full", "retmode": "xml", "ext": "xml"},
                {"rettype": "full", "retmode": "text", "ext": "txt"},
            ],
            "biosystems": [{"rettype": "xml", "retmode": "xml", "ext": "xml"}],
            "gds": [{"rettype": "summary", "retmode": "text", "ext": "txt"}],
            "gene": [
                {"rettype": "null", "retmode": "asn.1", "ext": "asn1"},
                {"rettype": "null", "retmode": "xml", "ext": "xml"},
                {"rettype": "gene_table", "retmode": "text", "ext": "gene_table"},
            ],
            "homologene": [
                {"rettype": "null", "retmode": "asn.1", "ext": "asn1"},
                {"rettype": "null", "retmode": "xml", "ext": "xml"},
                {
                    "rettype": "alignmentscores",
                    "retmode": "text",
                    "ext": "alignmentscores",
                },
                {"rettype": "fasta", "retmode": "text", "ext": "fasta"},
                {"rettype": "homologene", "retmode": "text", "ext": "homologene"},
            ],
            "mesh": [{"rettype": "full", "retmode": "text", "ext": "txt"}],
            "nlmcatalog": [
                {"rettype": "null", "retmode": "text", "ext": "txt"},
                {"rettype": "null", "retmode": "xml", "ext": "xml"},
            ],
            "nuccore": [
                {"rettype": "null", "retmode": "text", "ext": "txt"},
                {"rettype": "null", "retmode": "asn.1", "ext": "asn1"},
                {"rettype": "native", "retmode": "xml", "ext": "xml"},
                {"rettype": "acc", "retmode": "text", "ext": "acc"},
                {"rettype": "fasta", "retmode": "text", "ext": "fasta"},
                {"rettype": "fasta", "retmode": "xml", "ext": "fasta.xml"},
                {"rettype": "seqid", "retmode": "text", "ext": "seqid"},
                {"rettype": "gb", "retmode": "text", "ext": "gb"},
                {"rettype": "gb", "retmode": "xml", "ext": "gb.xml"},
                {"rettype": "gbc", "retmode": "xml", "ext": "gbc"},
                {"rettype": "ft", "retmode": "text", "ext": "ft"},
                {"rettype": "gbwithparts", "retmode": "text", "ext": "gbwithparts"},
                {"rettype": "fasta_cds_na", "retmode": "text", "ext": "fasta_cds_na"},
                {"rettype": "fasta_cds_aa", "retmode": "text", "ext": "fasta_cds_aa"},
            ],
            "nucest": [
                {"rettype": "null", "retmode": "text", "ext": "txt"},
                {"rettype": "null", "retmode": "asn.1", "ext": "asn1"},
                {"rettype": "native", "retmode": "xml", "ext": "xml"},
                {"rettype": "acc", "retmode": "text", "ext": "acc"},
                {"rettype": "fasta", "retmode": "text", "ext": "fasta"},
                {"rettype": "fasta", "retmode": "xml", "ext": "fasta.xml"},
                {"rettype": "seqid", "retmode": "text", "ext": "seqid"},
                {"rettype": "gb", "retmode": "text", "ext": "gb"},
                {"rettype": "gb", "retmode": "xml", "ext": "gb.xml"},
                {"rettype": "gbc", "retmode": "xml", "ext": "gbc"},
                {"rettype": "est", "retmode": "text", "ext": "est"},
            ],
            "nucgss": [
                {"rettype": "null", "retmode": "text", "ext": "txt"},
                {"rettype": "null", "retmode": "asn.1", "ext": "asn1"},
                {"rettype": "native", "retmode": "xml", "ext": "xml"},
                {"rettype": "acc", "retmode": "text", "ext": "acc"},
                {"rettype": "fasta", "retmode": "text", "ext": "fasta"},
                {"rettype": "fasta", "retmode": "xml", "ext": "fasta.xml"},
                {"rettype": "seqid", "retmode": "text", "ext": "seqid"},
                {"rettype": "gb", "retmode": "text", "ext": "gb"},
                {"rettype": "gb", "retmode": "xml", "ext": "gb.xml"},
                {"rettype": "gbc", "retmode": "xml", "ext": "gbc"},
                {"rettype": "gss", "retmode": "text", "ext": "gss"},
            ],
            "protein": [
                {"rettype": "null", "retmode": "text", "ext": "txt"},
                {"rettype": "null", "retmode": "asn.1", "ext": "asn1"},
                {"rettype": "native", "retmode": "xml", "ext": "xml"},
                {"rettype": "acc", "retmode": "text", "ext": "acc"},
                {"rettype": "fasta", "retmode": "text", "ext": "fasta"},
                {"rettype": "fasta", "retmode": "xml", "ext": "fasta.xml"},
                {"rettype": "seqid", "retmode": "text", "ext": "seqid"},
                {"rettype": "ft", "retmode": "text", "ext": "ft"},
                {"rettype": "gp", "retmode": "text", "ext": "gp"},
                {"rettype": "gp", "retmode": "xml", "ext": "gp.xml"},
                {"rettype": "gpc", "retmode": "xml", "ext": "gpc"},
                {"rettype": "ipg", "retmode": "xml", "ext": "xml"},
            ],
            "popset": [
                {"rettype": "null", "retmode": "text", "ext": "txt"},
                {"rettype": "null", "retmode": "asn.1", "ext": "asn1"},
                {"rettype": "native", "retmode": "xml", "ext": "xml"},
                {"rettype": "acc", "retmode": "text", "ext": "acc"},
                {"rettype": "fasta", "retmode": "text", "ext": "fasta"},
                {"rettype": "fasta", "retmode": "xml", "ext": "fasta.xml"},
                {"rettype": "seqid", "retmode": "text", "ext": "seqid"},
                {"rettype": "gb", "retmode": "text", "ext": "gb"},
                {"rettype": "gb", "retmode": "xml", "ext": "gb.xml"},
                {"rettype": "gbc", "retmode": "xml", "ext": "gbc"},
            ],
            "pmc": [
                {"rettype": "null", "retmode": "xml", "ext": "xml"},
                {"rettype": "medline", "retmode": "text", "ext": "medline"},
            ],
            "pubmed": [
                {"rettype": "null", "retmode": "asn.1", "ext": "asn1"},
                {"rettype": "null", "retmode": "xml", "ext": "xml"},
                {"rettype": "medline", "retmode": "text", "ext": "medline"},
                {"rettype": "uilist", "retmode": "text", "ext": "uilist"},
                {"rettype": "abstract", "retmode": "text", "ext": "abstract"},
            ],
            "sequences": [
                {"rettype": "null", "retmode": "text", "ext": "txt"},
                {"rettype": "acc", "retmode": "text", "ext": "acc"},
                {"rettype": "fasta", "retmode": "text", "ext": "fasta"},
                {"rettype": "seqid", "retmode": "text", "ext": "seqid"},
            ],
            "snp": [
                {"rettype": "null", "retmode": "asn.1", "ext": "asn1"},
                {"rettype": "null", "retmode": "xml", "ext": "xml"},
                {"rettype": "flt", "retmode": "text", "ext": "flt"},
                {"rettype": "fasta", "retmode": "text", "ext": "fasta"},
                {"rettype": "rsr", "retmode": "text", "ext": "rsr"},
                {"rettype": "ssexemplar", "retmode": "text", "ext": "ssexemplar"},
                {"rettype": "chr", "retmode": "text", "ext": "chr"},
                {"rettype": "docset", "retmode": "text", "ext": "docset"},
                {"rettype": "uilist", "retmode": "text", "ext": "uilist"},
                {"rettype": "uilist", "retmode": "xml", "ext": "uilist.xml"},
            ],
            "sra": [{"rettype": "full", "retmode": "xml", "ext": "xml"}],
            "taxonomy": [
                {"rettype": "null", "retmode": "xml", "ext": "xml"},
                {"rettype": "uilist", "retmode": "text", "ext": "uilist"},
                {"rettype": "uilist", "retmode": "xml", "ext": "uilist.xml"},
            ],
        }

    @property
    def valid_extensions(self):
        extensions = set()
        for db, db_options in self.efetch_options.items():
            for options in db_options:
                extensions |= set([options["ext"]])
        return list(extensions)

    def dbs_for_options(self, file_ext, rettype=None, retmode=None):
        possible_dbs = set()
        for db, db_options in self.efetch_options.items():
            for option_dict in db_options:
                if option_dict["ext"] == file_ext:
                    if retmode and option_dict["retmode"] != retmode:
                        continue
                    if rettype and option_dict["rettype"] != rettype:
                        continue
                    possible_dbs |= set([db])
                    break
        return possible_dbs

    def options_for_db_and_extension(self, db, file_ext, rettype=None, retmode=None):
        possible_options = []
        assert file_ext, "file_ext must be defined"

        if not self.is_valid_db(db):
            raise NCBIFileException(
                "DB specified is not valid. Options include: {dbs}".format(
                    dbs=", ".join(self.valid_dbs)
                )
            )

        db_options = self.efetch_options[db]
        for opt in db_options:
            if file_ext == opt["ext"]:
                if retmode and opt["retmode"] != retmode:
                    continue
                if rettype and opt["rettype"] != rettype:
                    continue
                possible_options.append(opt)

        return possible_options

    def guess_db_options_for_extension(
        self, file_ext, db=None, rettype=None, retmode=None
    ):
        if db and rettype and retmode:
            if self.is_valid_db_request(db, rettype, retmode):
                request_options = {}
                request_options["db"] = db
                request_options["rettype"] = rettype
                request_options["retmode"] = retmode
                request_options["ext"] = file_ext
                return request_options

        possible_dbs = [db] if db else self.dbs_for_options(file_ext, rettype, retmode)

        if len(possible_dbs) > 1:
            raise NCBIFileException(
                'Ambigious db for file extension specified: "{}"; possible databases include: {}'.format(
                    file_ext, ", ".join(list(possible_dbs))
                )
            )
        elif len(possible_dbs) == 1:
            likely_db = possible_dbs.pop()

            likely_options = self.options_for_db_and_extension(
                likely_db, file_ext, rettype, retmode
            )
            if len(likely_options) == 1:
                request_options = {}
                request_options["db"] = likely_db
                request_options["rettype"] = likely_options[0]["rettype"]
                request_options["retmode"] = likely_options[0]["retmode"]
                request_options["ext"] = likely_options[0]["ext"]
                return request_options
            elif len(likely_options) > 1:
                raise NCBIFileException(
                    "Please clarify the rettype and retmode. Multiple request types are possible for the file extension ({}) specified: {}".format(
                        file_ext, likely_options
                    )
                )
            else:
                raise NCBIFileException(
                    "No request options found. Please check the file extension ({}), db ({}), rettype ({}), and retmode ({}) specified.".format(
                        file_ext, db, rettype, retmode
                    )
                )

    def is_valid_db_request(self, db, rettype, retmode):
        if not self.is_valid_db(db):
            raise NCBIFileException(
                "DB specified is not valid. Options include: {dbs}".format(
                    dbs=", ".join(self.valid_dbs)
                )
            )
        db_options = self.efetch_options[db]
        for opt in db_options:
            if opt["rettype"] == rettype and opt["retmode"] == retmode:
                return True
        return False

    @property
    def valid_dbs(self):
        return self.efetch_options.keys()

    def is_valid_db(self, db):
        return db in self.valid_dbs

    def parse_accession_str(self, id_str):
        """
            This tries to match an NCBI accession as defined here:
                https://www.ncbi.nlm.nih.gov/Sequin/acc.html
        """
        m = re.search(
            r"(?P<accession>(?:[a-zA-Z]{1,6}|NW_|NC_|NM_|NR_)\d{1,10})(?:\.(?P<version>\d+))?(?:\.(?P<file_ext>\S+))?.*",
            id_str,
        )
        accession, version, file_ext = ("", "", "")
        if m:
            accession = m.group("accession")
            version = m.group("version")
            file_ext = m.group("file_ext")
        assert (
            file_ext
        ), "file_ext must be defined: {}.{}.<file_ext>. Possible values include: {}".format(
            accession, version, ", ".join(list(self.valid_extensions))
        )
        assert version, "version must be defined: {}.<version>.{}".format(
            accession, file_ext
        )

        return accession, version, file_ext

    @staticmethod
    def _seq_chunks(seq, n):
        # https://stackoverflow.com/a/312464/190597 (Ned Batchelder)
        """ Yield successive n-sized chunks from seq."""
        for i in range(0, len(seq), n):
            yield seq[i : i + n]

    def _esummary_and_parse(
        self,
        accession,
        xpath_selector,
        db="nuccore",
        return_type=int,
        raise_on_failure=True,
        retmode="xml",
        **kwargs
    ):
        result = self.entrez.esummary(db=db, id=accession, **kwargs)

        root = ET.fromstring(result.read())
        nodes = root.findall(xpath_selector)

        retval = 0
        if len(nodes):
            retval = return_type(nodes[0].text)
        else:
            if raise_on_failure:
                raise NCBIFileException("The esummary query failed.")

        return retval

    def exists(self, accession, db="nuccore"):
        result = self.entrez.esearch(db=db, term=accession, rettype="count")

        root = ET.fromstring(result.read())
        nodes = root.findall(".//Count")

        count = 0
        if len(nodes):
            count = int(nodes[0].text)
        else:
            raise NCBIFileException("The esummary query failed.")

        if count == 1:
            return True
        else:
            logger.warning(
                'The accession specified, "{acc}", could not be found in the database "{db}".\nConsider if you may need to specify a different database via "db=<db_id>".'.format(
                    acc=accession, db=db
                )
            )
            return False

    def size(self, accession, db="nuccore"):
        return self._esummary_and_parse(accession, ".//*[@Name='Length']", db=db)

    def mtime(self, accession, db="nuccore"):
        update_date = self._esummary_and_parse(
            accession, ".//Item[@Name='UpdateDate']", db=db, return_type=str
        )

        pattern = "%Y/%m/%d"
        epoch_update_date = int(time.mktime(time.strptime(update_date, pattern)))

        return epoch_update_date

    def fetch_from_ncbi(
        self,
        accession_list,
        destination_dir,
        force_overwrite=False,
        rettype="fasta",
        retmode="text",
        file_ext=None,
        combined_file_prefix=None,
        remove_separate_files=False,
        chunk_size=1,
        db="nuccore",
        **kwargs
    ):
        """
            This function downloads and saves files from NCBI.
            Adapted in part from the BSD-licensed code here:
              https://github.com/broadinstitute/viral-ngs/blob/master/util/genbank.py
        """

        max_chunk_size = 500

        # Conform to NCBI retreival guidelines by chunking into 500-accession chunks if
        # >500 accessions are specified and chunk_size is set to 1
        # Also clamp chunk size to 500 if the user specified a larger value.
        if chunk_size > max_chunk_size or (
            len(accession_list) > max_chunk_size and chunk_size == 1
        ):
            chunk_size = max_chunk_size

        outEx = {"fasta": "fasta", "ft": "tbl", "gb": "gbk"}

        output_directory = os.path.abspath(os.path.expanduser(destination_dir))

        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        output_extension = str(file_ext)

        # ensure the extension starts with a ".", also allowing for passed-in
        # extensions that already have it
        if output_extension[:1] != ".":
            output_extension = "." + output_extension

        logger.info(
            "Fetching {} entries from NCBI: {}\n".format(
                str(len(accession_list)), ", ".join(accession_list[:10])
            )
        )
        output_files = []

        for chunk_num, chunk in enumerate(self._seq_chunks(accession_list, chunk_size)):
            # sleep to throttle requests to 2 per second per NCBI guidelines:
            #   https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
            time.sleep(0.5)
            acc_string = ",".join(chunk)

            # if the filename would be longer than Linux allows, simply say "chunk-chunk_num"
            if len(acc_string) + len(output_extension) <= 254:
                output_file_path = os.path.join(
                    output_directory, acc_string + output_extension
                )
            else:
                output_file_path = os.path.join(
                    output_directory, "chunk-{}".format(chunk_num) + output_extension
                )

            if not force_overwrite:
                logger.info("not overwriting, checking for existence")
                assert not os.path.exists(output_file_path), (
                    """File %s already exists. Consider removing
                    this file or specifying a different output directory. The files for the accessions specified
                    can be overwritten if you add force_overwrite flag. Processing aborted."""
                    % output_file_path
                )

            try_count = 1
            while True:
                try:
                    logger.info(
                        "Fetching file {}: {}, try #{}".format(
                            chunk_num + 1, acc_string, try_count
                        )
                    )
                    handle = self.entrez.efetch(
                        db=db, rettype=rettype, retmode=retmode, id=acc_string, **kwargs
                    )

                    with open(output_file_path, "w") as outf:
                        for line in handle:
                            outf.write(line)
                    output_files.append(output_file_path)
                except IOError:

                    logger.warning(
                        "Error fetching file {}: {}, try #{} probably because NCBI is too busy.".format(
                            chunk_num + 1, acc_string, try_count
                        )
                    )

                    try_count += 1
                    if try_count > 4:
                        logger.warning("Tried too many times. Aborting.")
                        raise

                    # if the fetch failed, wait a few seconds and try again.
                    logger.info("Waiting and retrying...")
                    time.sleep(2)

                    continue
                break

        # assert that we are not trying to remove the intermediate files without writing a combined file
        if remove_separate_files:
            assert combined_file_prefix, """The intermediate files
                can only be removed if a combined file is written via combined_file_prefix"""

        # build a path to the combined genome file
        if combined_file_prefix:
            concatenated_genome_file_path = os.path.join(
                output_directory, combined_file_prefix + output_extension
            )

            if not force_overwrite:
                assert not os.path.exists(concatenated_genome_file_path), (
                    """File %s already exists. Consider removing
                    this file or specifying a different output directory. The files for the accessions specified
                    can be overwritten if you add force_overwrite flag. Processing aborted."""
                    % output_file_path
                )

            # concatenate the files together into one genome file
            with open(concatenated_genome_file_path, "w") as outfile:
                for file_path in output_files:
                    with open(file_path) as infile:
                        for line in infile:
                            outfile.write(line)

            # if the option is specified, remove the intermediate fasta files
            if remove_separate_files:
                while len(output_files) > 0:
                    os.unlink(output_files.pop())

            # add the combined file to the list of files returned
            output_files.append(concatenated_genome_file_path)

        # return list of files
        return output_files

    def search(self, query, *args, db="nuccore", idtype="acc", **kwargs):
        # enforce JSON return mode
        kwargs["retmode"] = "json"

        # if the user specifies retmax, use it and limit there
        # otherwise page 200 at a time and return all
        if "retmax" not in kwargs or not kwargs["retmax"]:
            kwargs["retmax"] = 100000
            return_all = True
        else:
            return_all = False

        kwargs["retstart"] = kwargs.get("retstart", 0)

        def esearch_json(term, *args, **kwargs):
            handle = self.entrez.esearch(term=term, *args, **kwargs)
            json_result = json.loads(handle.read())
            return json_result

        def result_ids(json):
            if (
                "esearchresult" in json_results
                and "idlist" in json_results["esearchresult"]
            ):
                return json_results["esearchresult"]["idlist"]
            else:
                raise NCBIFileException("ESearch error")

        has_more = True

        while has_more:
            json_results = esearch_json(
                term=query, *args, db=db, idtype=idtype, **kwargs
            )

            for acc in result_ids(json_results):
                yield acc

            if return_all and (
                "count" in json_results["esearchresult"]
                and int(json_results["esearchresult"]["count"])
                > kwargs["retmax"] + kwargs["retstart"]
            ):
                kwargs["retstart"] += kwargs["retmax"]
                # sleep to throttle requests to <2 per second per NCBI guidelines:
                #   https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
                time.sleep(0.5)
            else:
                has_more = False
