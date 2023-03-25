__author__ = "Simone Baffelli"
__copyright__ = "Copyright 2022, Simone Baffelli"
__email__ = "simone.baffelli@gmail.com"
__license__ = "MIT"

import os
import re
import time
from pathlib import Path
from re import Match
from typing import Optional, Tuple

import sqlalchemy.databases
from sqlalchemy import MetaData, Table, create_engine, engine, sql
from sqlalchemy.engine.url import URL

from snakemake.exceptions import (HTTPFileException, RuleException,
                                  WorkflowError)
from snakemake.logging import logger
# module-specific
from snakemake.remote import (AbstractRemoteObject, AbstractRemoteProvider,
                              DomainObject)

# def row_to_dict(cur):
# return dict(zip([t[0] for t in row.cursor_description], row))


class SQLFileException(RuleException):
    def __init__(self, msg, lineno=None, snakefile=None):
        super().__init__(msg, lineno=lineno, snakefile=snakefile)


# class MySQLConnectionManager():
#     def __init__(self, *args, **kwargs):
#         self.args = args
#         self.kwargs = kwargs
#         self.connection = None
#         self.eng = None

#     def create_conn(self):
#             return mysql.connector.connect(*self.args, **self.kwargs)

#     def __enter__(self):
#         self.eng  = create_engine("mysql+mysqlconnector://", creator=self.create_conn)
#         self.connection = self.eng.connect()
#         return self.connection

#     def __exit__(self, exc_type, exc_value, exc_traceback):
#         self.connection.close()

#     def cursor(self):
#         return ManagedCursor(self.connection)


# class ManagedCursor():
#     def __init__(self, conn, **kwargs):
#         self.conn = conn
#         self.cur = None

#     def __enter__(self):
#         self.cur = self.conn.cursor(dictionary=True)
#         return self

#     def __exit__(self, exc_type, exc_value, exc_traceback):
#         if not self.conn.is_connected():
#             self.conn.reconnect()
#         self.cur.close()

#     def cursor(self):
#         if self.cur:
#             return self.cur
#         else:
#             return self.conn.cursor(dictionary=True)

#     def close(self):
#         if self.cur:
#             self.cur.close()

#     def execute(self, *args, **kwargs):
#         if self.cur:
#             try:
#                 self.cur.execute(*args, **kwargs)
#             except msqe.ProgrammingError as e:
#                 print(e)
#                 raise e
#         else:
#             raise(ConnectionError('There is no open MySql Connection'))


#     def fetchone(self):
#         return self.cur.fetchone()


# def backup_database(**kwargs, path):
#     cmd= sql.SQL("""\
#     USE {db};
#     BACKUP DATABASE {db}
#     TO DISK = 'c:\tmp\{db}.bak'
#     WITH FORMAT,
#         MEDIANAME = 'SQLServerBackups',
#         NAME = 'Full Backup of {db}';
#     GO
#     """.fomat(db=sql.identifier(kwargs.get(database))))
#     return cmd


def run_query_on_connection(query: str, eng: engine):
    try:
        with eng.connect() as conn:
            conn.execute(query)
    except Exception as e:
        print(e)
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


# class BulkLoggeddDbSession:
#     """
#     Class to represent an bulk logged db session
#     This is a context manager which switches the
#     logging behavior of the db in the context
#     """
#     def __init__(self, *args, **kwargs):
#         self.params = kwargs
#         self._conn = mysql.connector.connect(*args, as_dict=True,**kwargs)

#     def __enter__(self):
#         alter_query =\
#         f"""
#         ALTER DATABASE [{self.params['database']}] SET RECOVERY BULK_LOGGED;
#         """
#         with self._conn.cursor() as cur:
#             cur.execute(alter_query)

#     def __exit__(self, type, value, traceback):
#         alter_query =\
#         f"""
#         ALTER DATABASE [{self.params['database']}] SET RECOVERY FULL;
#         """
#         with self._conn.cursor() as cur:
#             cur.execute(alter_query)

#     def run_query(self, query):
#         run_query_on_connection(query, self._conn)


# class ExplicitelyCommitedTransaction:
#     """
#     Implement a context manager in which SQL transactions are
#     only commited after leaving the context
#     """
#     def __init__(self, *args, **kwargs):
#         self.params = kwargs
#         #self._conn = pymssql.connect(*args, as_dict=True,**kwargs)


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
        return "sql://"

    @property
    def available_protocols(self):
        """List of valid protocols for this remote provider."""
        db_names = sqlalchemy.databases.__all__
        return ["sql://", "postgresql://", "mysql://", "sqlite://"]

    def connect(self, **kwargs):
        print("Connecting")
        return self._sqlc.connect(**kwargs)


class RemoteObject(AbstractRemoteObject):
    """
    This is a class to interact with an SQL server as if it were a file
    """

    def __init__(
        self,
        *args,
        keep_local:bool=False,
        stay_on_remote:bool=True,
        provider=None,
        ancient:bool=False,
        time_query:str=None
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
        md = MetaData(bind=self._sqlc)
        md.reflect()
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
        return f"{self.schema}.{self.table}"

    def check_existence(self) -> bool:
        """
        Check if the table exists by querying the information schema
        """
        breakpoint()
        return self._sqlc.has_table(self.table, schema=self.schema)

    def exists(self):
        """
        This is used by Snakemake to check if the "file" represented by the table exists

        """
        return self.check_existence()

    def mtime(self):
        if self.exists() and not self.ancient:
            return self.last_modified()
        elif self.exists() and self.ancient:
            return self.last_modified() - 1000000000
        else:
            raise SQLFileException(
                f"The file cannot be found in {self.fully_qualified_table}"
            )

    def is_newer(self, time):
        if not self.ancient:
            return self.mtime() > time
        else:
            return False
        
    def get_table(self) -> Table:
        return self._md.tables[self.table]

    def modified_query(self):
        tb = self.get_table()
        tb.select(sql.expression.text(self.time_query).)

    def last_modified(self):
        breakpoint()
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
        size_query = f"""
        SELECT 
            COUNT(*) AS n 
        FROM {self.fully_qualified_table} 
        """.strip()
        if self.check_existence():
            with self._sqlc.connect() as cur:
                res = cur.execute(size_query)
                sz = res.fetchone()[0]
                res.fetchall()
            return sz
        else:
            raise SQLFileException(
                f"The file cannot be found in {self.fully_qualified_table}"
            )

    def remove(self):
        pass
        rm_query = f"""
        DROP TABLE IF EXISTS {self.fully_qualified_table};
        """.strip()
        if self.check_existence():
            try:
                with self._sqlc.connection() as cur:
                    res = cur.execute(rm_query)
                    res.fetchall()
            except:
                return True
        return True

    def download(self, make_dest_dirs=True):
        if self.check_existence():
            if make_dest_dirs:
                os.makedirs(os.path.dirname(self.local_path), exist_ok=True)
                Path(self.local_path).touch()
                os.sync()  # ensure flush to disk
            else:
                raise SQLFileException(
                    f"The file cannot be found in {self.fully_qualified_table}"
                )
