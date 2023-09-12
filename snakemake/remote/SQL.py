__author__ = "Simone Baffelli"
__copyright__ = "Copyright 2022, Simone Baffelli"
__email__ = "simone.baffelli@gmail.com"
__license__ = "MIT"

import os
import re
import time
from pathlib import Path
from re import Match
from typing import Optional, Tuple, Callable

import sqlalchemy.dialects as databases
from sqlalchemy import MetaData, Table, create_engine, engine, sql, inspect
from sqlalchemy.engine import Inspector
from sqlalchemy.engine.url import URL

from snakemake.exceptions import HTTPFileException, RuleException, WorkflowError
from snakemake.logging import logger

# module-specific
from snakemake.remote import AbstractRemoteObject, AbstractRemoteProvider, DomainObject


import datetime as dt


class SQLFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)


def run_query_on_connection(query: str, eng: engine):
    try:
        with eng.connect() as conn:
            conn.execute(query)
    except Exception as e:
        raise (e)


def run_query_from_file_on_connection(query, conn, multi=True):
    with open(query, "r") as inf:
        fc = inf.readlines()
        full_query = "".join(fc)
        queries = [q for q in full_query.split(";") if q != ""]
        res = [conn.execute(q) for q in queries if q != "\n"]
        valids = [r.fetchall() for r in res if r.rowcount > 0 and r.returns_rows]
        [r.close() for r in valids]
        return valids


class RemoteProvider(AbstractRemoteProvider):
    def __init__(
        self, *args, keep_local=False, stay_on_remote=True, is_default=False, **kwargs
    ):
        super(RemoteProvider, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            is_default=is_default,
            **kwargs,
        )
        # Create sqlalchemy engine
        self._sqlc = engine.create_engine(*args, **kwargs)

    @property
    def remote_interface(self):
        return self._sqlc

    @property
    def default_protocol(self):
        """The protocol that is prepended to the path when no protocol is specified."""
        return "jdbc://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        db_names = ["{name}://" for name in databases.__all__]
        return db_names + ["sqlite://", "duckdb://"]

    def connect(self, **kwargs):
        return self._sqlc.connect(**kwargs)


class RemoteObject(AbstractRemoteObject):
    """
    This is a class to interact with an SQL server as if it were a file
    """

    _sqlc: engine.Engine

    def __init__(
        self,
        *args,
        keep_local: bool = False,
        stay_on_remote: bool = True,
        provider=None,
        ancient: bool = False,
        time_query: Optional[str] = None,
        **kwargs,
    ):
        super(RemoteObject, self).__init__(
            *args,
            keep_local=keep_local,
            stay_on_remote=stay_on_remote,
            provider=provider,
            **kwargs,
        )
        if provider:
            self._sqlc = provider.remote_interface
        else:
            self._sqlc = args[0]
        md = MetaData()
        md.reflect(bind=self._sqlc)
        self._md = md
        self.ancient = ancient
        self.time_query = time_query

    @property
    def _matched_schema_table(self) -> Tuple[Optional[str], Optional[str]]:
        matches = re.search(
            r"(((\[?)(?P<schema>\w+)(\]?))\.)?(((\[?)(?P<table>\w+)(\]?)))",
            self.local_file(),
        )
        if matches is not None:
            schema = matches.group("schema")
            table = matches.group("table")
            return (schema, table)
        else:
            return (None, None)

    @property
    def schema(self) -> Optional[str]:
        schema, table = self._matched_schema_table
        if table is not None:
            return schema

    @property
    def table(self) -> str:
        schema, table = self._matched_schema_table
        if table is not None:
            return table
        else:
            return schema

    @property
    def fully_qualified_table(self) -> str:
        if self.schema:
            return f"{self.schema}.{self.table}"
        else:
            return self.table

    def check_existence(self) -> bool:
        """
        Check if the table exists by querying the information schema
        """
        insp = inspect(self._sqlc)
        return insp.has_table(self.table, schema=self.schema)

    def exists(self) -> bool:
        """

        This is used by Snakemake to check if the "file" represented by the table exists

        """
        return self.check_existence()

    def mtime(self):
        if self.exists() and not self.ancient:
            return self.last_modified()
        elif self.exists() and self.ancient:
            return float("-Inf")
        else:
            return float("-Inf")

    def is_newer(self, time: float) -> bool:
        if not self.ancient:
            return self.mtime() > time
        else:
            return False

    def get_table(self) -> Table:
        self._md.reflect(bind=self._sqlc)
        return self._md.tables[self.table]

    def modified_query(self) -> float:
        tb = self.get_table()
        with self._sqlc.connect() as cur:
            qr = sql.select(sql.literal_column((self.time_query))).select_from(tb)
            result = cur.execute(qr)
        dates = [r for r, *rest in result.all()]
        if len(dates) == 1 and dates[0]:
            return float(dates[0])
        elif len(dates) > 1:
            raise ValueError("The query returned more than one row")
        else:
            return float("-Inf")

    def last_modified(self) -> float:
        if self.time_query:
            return self.modified_query()
        else:
            mod_query = f"""
            SELECT
                TABLE_NAME,
                TABLE_SCHEMA,
                COALESCE(UPDATE_TIME, CREATE_TIME)  AS modification_date
            FROM INFORMATION_SCHEMA.TABLES 
            WHERE TABLE_TYPE IN ('BASE TABLE','VIEW') AND TABLE_NAME = \'{self.table}\' AND TABLE_SCHEMA = \'{self.schema}\'
            """.strip()
            with self._sqlc.connect() as cur:
                res = cur.execute(mod_query)
                rw = res.fetchall()
            if len(rw) > 0:
                md = rw[0]["modification_date"]
                if md:
                    return md.timestamp()

    def size(self):
        tb = self.get_table()
        if self.check_existence():
            with self._sqlc.connect() as cur:
                res = cur.execute(sql.select(sql.func.count()).select_from(tb))
                sz = res.fetchone()[0]
                res.fetchall()
            return sz
        else:
            raise SQLFileException(
                f"The file cannot be found in {self.fully_qualified_table}"
            )

    def remove(self):
        if self.check_existence():
            try:
                with self._sqlc.connect() as cur:
                    tb = self.get_table()
                    tb.drop(bind=cur)
            except:
                raise SQLFileException(
                    f"Cannot drop table {self.fully_qualified_table}"
                )
        return True

    def _download(self, make_dest_dirs=True):
        if self.check_existence():
            if make_dest_dirs:
                os.makedirs(os.path.dirname(self.local_path), exist_ok=True)
                Path(self.local_path).touch()
                os.sync()  # ensure flush to disk
            else:
                raise SQLFileException(
                    f"The file cannot be found in {self.fully_qualified_table}"
                )
