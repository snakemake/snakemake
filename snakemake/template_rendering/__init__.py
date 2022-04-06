from abc import ABC, abstractmethod

from snakemake.exceptions import WorkflowError


class TemplateRenderer(ABC):
    def __init__(self, input, output, params, wildcards, config):
        if len(output) != 1:
            raise ValueError(
                "More than one output file specified for template_engine rule."
            )
        if len(input) != 1:
            if "template" not in input.keys():
                raise ValueError(
                    "More than one input file specified for template engine rule, but no "
                    "input file named as 'template'."
                )
            else:
                self.input_file = input.template
        else:
            self.input_file = input[0]

        self.input = input
        self.output_file = output[0]
        self.params = params
        self.wildcards = wildcards
        self.config = config

    @property
    def variables(self):
        return {
            "params": self.params,
            "wildcards": self.wildcards,
            "config": self.config,
            "input": self.input,
        }

    @abstractmethod
    def render(self):
        ...


def render_template(engine, input, output, params, wildcards, config, rule):
    try:
        if engine == "yte":
            from snakemake.template_rendering.yte import YteRenderer

            return YteRenderer(input, output, params, wildcards, config).render()
        elif engine == "jinja2":
            from snakemake.template_rendering.jinja2 import Jinja2Renderer

            return Jinja2Renderer(input, output, params, wildcards, config).render()
        else:
            raise WorkflowError(
                f"Unsupported template engine {engine} in rule {rule}. "
                "So far, only yte and jinja2 are supported."
            )
    except Exception as e:
        raise WorkflowError(f"Error rendering template in rule {rule}.", e)
