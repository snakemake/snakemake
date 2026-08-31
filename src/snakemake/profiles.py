from collections import OrderedDict
import os
from pathlib import Path
from configargparse import YAMLConfigFileParser, ConfigFileParserException

from snakemake_interface_common.exceptions import WorkflowError


class ProfileConfigFileParser(YAMLConfigFileParser):
    def parse(self, stream):
        # taken from configargparse and modified to add special handling for key-value pairs
        import yaml

        profile_name = getattr(stream, "name", "<stream>")
        try:
            text = stream.read()
            parsed_obj = yaml.load(text, Loader=yaml.FullLoader)

            if parsed_obj.pop("__use_yte__", False):
                try:
                    import yte
                except ImportError as e:
                    raise WorkflowError(
                        f"Config file '{profile_name}' requires the missing package 'yte'. "
                        f"Please install it with 'pip install yte'."
                    ) from e

                parsed_obj = yte.process_yaml(text, require_use_yte=True)

            if not isinstance(parsed_obj, dict):
                raise KeyError(
                    f"It doesn't appear to contain 'key: value' pairs (aka. a YAML mapping). "
                    f"Got '{type(parsed_obj).__name__}' instead of 'dict'."
                )
        except WorkflowError:
            raise
        except Exception as e:
            raise ConfigFileParserException(
                f"Couldn't parse config file '{profile_name}'."
            ) from e

        def format_val(val):
            def repr_if_numeric(val, numtype):
                try:
                    numtype(val)
                    return repr(val)
                except ValueError:
                    return None

            if isinstance(val, str):
                # when the value is an explicit string enclosing int or float, we have
                # to quote it to prevent an evaluation as integer in the next pass
                return repr_if_numeric(val, int) or repr_if_numeric(val, float) or val
            return val

        def format_one_level_dict(d):
            return [f"{key}={format_val(val)}" for key, val in d.items()]

        def format_two_level_dict(d, item: str):
            if not all(isinstance(val, dict) for val in d.values()):
                raise WorkflowError(
                    f"Invalid {item} format in profile. Expected two-level mapping, got {d}"
                )
            return [
                f"{key}:{key2}={format_val(val2)}"
                for key, val in d.items()
                for key2, val2 in val.items()
            ]

        def resolve_path(raw):
            """
            Adjust path if it exists in the profile dir.
            Otherwise value is not a file or not existing in the profile dir.

            Cache files in the profile directory for performance.
            Otherwise fall back to search in the filesystem.
            """
            value = os.path.expanduser(os.path.expandvars(str(raw)))
            file = profile_dir / value
            if file.parent == profile_dir:
                if not profile_cache:
                    profile_cache.update(os.listdir(profile_dir))
                if value in profile_cache:
                    return str(file)
            elif file.exists():
                return str(file)

        profile_dir = Path(profile_name).parent
        profile_cache = set()
        result = OrderedDict()
        for key, value in parsed_obj.items():
            if isinstance(value, list):
                result[key] = value
            elif value is None:
                continue
            else:
                # special handling for simplified pure YAML syntax for key-value CLI arguments like --resources
                if isinstance(value, (dict, OrderedDict)):
                    if key in (
                        "set-threads",
                        "resources",
                        "set-scatter",
                        "batch",
                        "set-resource-scopes",
                        "default-resources",
                        "config",
                        "groups",
                        "group-components",
                        "consider-ancient",
                    ):
                        result[key] = format_one_level_dict(value)
                    elif key == "set-resources":
                        result[key] = format_two_level_dict(value, "set-resources")
                else:
                    value = resolve_path(value)
                    if value:
                        result[key] = value

        return result
