from abc import ABC, abstractmethod
from collections import namedtuple
from collections.abc import Mapping, Callable
from functools import partial
import inspect
import os
import re
from typing import List, Optional, Union

from snakemake.common import async_run
import snakemake.io
import snakemake.utils
from snakemake.exceptions import LookupError
from snakemake_interface_common.exceptions import WorkflowError


# Alias for expand that provides a more intuitive name for the use case of
# collecting files from previous jobs.
collect = snakemake.io.expand
