from typing import List
from snakemake_interface_common.settings import SettingsEnumBase


class RerunTrigger(SettingsEnumBase):
    MTIME = 0
    PARAMS = 1
    INPUT = 2
    SOFTWARE_ENV = 3
    CODE = 4


class ChangeType(SettingsEnumBase):
    CODE = 0
    INPUT = 1
    PARAMS = 2


class CondaCleanupPkgs(SettingsEnumBase):
    TARBALLS = 0
    CACHE = 1


class Quietness(SettingsEnumBase):
    RULES = 0
    PROGRESS = 1
    ALL = 2
    HOST = 3
