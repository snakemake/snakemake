import argparse
import collections
import dataclasses

import configargparse


class ArgumentParser(configargparse.ArgumentParser):
    def add_argument(
        self,
        *args,
        parse_func=None,
        **kwargs,
    ):
        if parse_func is not None:
            register_parser_action(parse_func, kwargs)
        super().add_argument(*args, **kwargs)

    def add_argument_group(self, *args, **kwargs):
        group = ArgumentGroup(self, *args, **kwargs)
        self._action_groups.append(group)
        return group


class ArgumentGroup(argparse._ArgumentGroup):
    def add_argument(
        self,
        *args,
        parse_func=None,
        **kwargs,
    ):
        if parse_func is not None:
            register_parser_action(parse_func, kwargs)
        super().add_argument(*args, **kwargs)


def register_parser_action(parse_func, kwargs):
    if "action" in kwargs:
        raise ValueError(
            "Cannot specify action if parser argument is provided to add_argument."
        )

    class ParserAction(argparse._StoreAction):
        def __init__(self, *args, **kwargs):
            if "parser" in kwargs:
                del kwargs["parse_func"]
            super().__init__(*args, **kwargs)

        def __call__(
            self,
            parser,
            namespace,
            values,
            option_string=None,
        ):
            parsed = parse_func(values)
            setattr(namespace, self.dest, parsed)

    kwargs["action"] = ParserAction


class ArgumentDefaultsHelpFormatter(argparse.HelpFormatter):
    """Help message formatter which adds default values to argument help.

    Like argparse.ArgumentDefaultsHelpFormatter, but doesn't print
    None/dataclasses._MISSING_TYPE/etc.
    """

    def _get_help_string(self, action):
        if (
            (
                action.option_strings
                or action.nargs in [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
            )
            and action.default not in (None, "", set(), argparse.SUPPRESS)
            and not isinstance(action.default, dataclasses._MISSING_TYPE)
        ):
            if isinstance(action.default, collections.abc.Iterable) and not isinstance(
                action.default, str
            ):

                if isinstance(action.default, (frozenset, set)):
                    default = sorted(map(str, action.default))
                else:
                    default = map(str, action.default)
                return action.help + f" (default: {' '.join(default)})"
            else:
                return action.help + " (default: %(default)s)"
        else:
            return action.help
