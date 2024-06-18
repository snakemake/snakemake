import argparse

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
