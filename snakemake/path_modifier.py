__authors__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
from snakemake.exceptions import WorkflowError
from snakemake.io import is_flagged, AnnotatedString, flag, get_flag_value


PATH_MODIFIER_FLAG = "path_modified"


class PathModifier:
    def __init__(self, replace_prefix: dict, workflow):
        self.skip_properties = set()
        self.workflow = workflow

        self.trie = None
        if replace_prefix:
            import datrie

            self.trie = datrie.Trie(
                "".join(set(char for prefix in replace_prefix for char in prefix))
            )
            for prefix, replacement in replace_prefix.items():
                self.trie[prefix] = replacement

    def modify(self, path, property=None):
        if get_flag_value(path, PATH_MODIFIER_FLAG) is self:
            # Path has been modified before and is reused now, no need to modify again.
            return path

        modified_path = self.apply_default_remote(self.replace_prefix(path, property))
        if modified_path == path:
            # nothing has changed
            return path

        # Important, update with previous flags in case of AnnotatedString #596
        if hasattr(path, "flags"):
            if not hasattr(modified_path, "flags"):
                modified_path = AnnotatedString(modified_path)
            modified_path.flags.update(path.flags)
        # Flag the path as modified and return.
        modified_path = flag(modified_path, PATH_MODIFIER_FLAG, self)
        return modified_path

    def replace_prefix(self, path, property=None):
        if self.trie is None or property in self.skip_properties:
            # no replacement
            return path
        prefixes = self.trie.prefix_items(str(path))
        if len(prefixes) > 1:
            # ambiguous prefixes
            raise WorkflowError(
                "Multiple prefixes ({}) match the path {}. Make sure that the replace_prefix statement "
                "in your module definition does not yield ambiguous matches.".format(
                    ", ".join(prefix[0] for prefix in prefixes), path
                )
            )
        elif prefixes:
            # replace prefix
            prefix, replacement = prefixes[0]
            return replacement + path[len(prefix) :]
        else:
            # no matching prefix
            return path

    def apply_default_remote(self, path):
        """Apply the defined default remote provider to the given path and return the updated _IOFile.
        Asserts that default remote provider is defined.
        """

        def is_annotated_callable(value):
            if isinstance(value, AnnotatedString):
                return bool(value.callable)

        if (
            self.workflow.default_remote_provider is None
            or is_flagged(path, "remote_object")
            or is_flagged(path, "local")
            or is_annotated_callable(path)
        ):
            # no default remote needed
            return path

        # This will convert any AnnotatedString to str
        fullpath = "{}/{}".format(self.workflow.default_remote_prefix, path)
        fullpath = os.path.normpath(fullpath)
        remote = self.workflow.default_remote_provider.remote(fullpath)
        return remote

    @property
    def modifies_prefixes(self):
        return self.trie is not None
