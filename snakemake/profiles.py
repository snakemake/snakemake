from collections import OrderedDict
from configargparse import YAMLConfigFileParser, ConfigFileParserException


class ProfileConfigFileParser(YAMLConfigFileParser):
    def parse(self, stream):
        # taken from configargparse and modified to add special handling for key-value pairs
        import yte

        try:
            parsed_obj = yte.process_yaml(stream, require_use_yte=True)
        except Exception as e:
            raise ConfigFileParserException("Couldn't parse config file: %s" % e)

        if not isinstance(parsed_obj, dict):
            raise ConfigFileParserException(
                "The config file doesn't appear to "
                "contain 'key: value' pairs (aka. a YAML mapping). "
                "yaml.load('%s') returned type '%s' instead of 'dict'."
                % (getattr(stream, "name", "stream"), type(parsed_obj).__name__)
            )

        def format_one_level_dict(d):
            return [f"{key}={val}" for key, val in d.items()]

        def format_two_level_dict(d):
            return [
                f"{key}:{key2}={val2}"
                for key, val in d.items()
                for key2, val2 in val.items()
            ]

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
                        "wms-monitor-arg",
                    ):
                        result[key] = format_one_level_dict(value)
                    if key in ("set-resources"):
                        result[key] = format_two_level_dict(value)
                else:
                    result[key] = str(value)

        return result
