from abc import ABC, abstractmethod


class TemplateRenderer(ABC):
    def __init__(self, input, output, params, wildcards, config):
        assert len(input) == 1
        assert len(output) == 1

        self.input_file = input[0]
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
        }

    @abstractmethod
    def render(self):
        ...


def render_template(engine, input, output, params, wildcards, config):
    if engine == "yte":
        from snakemake.template_rendering.yte import YteRenderer

        return YteRenderer(input, output, params, wildcards, config).render()
    elif engine == "jinja2":
        from snakemake.template_rendering.jinja2 import Jinja2Renderer

        return Jinja2Renderer(input, output, params, wildcards, config).render()
    else:
        raise WorkflowError(
            f"Unsupported template engine {engine}. "
            "So far, only yte and jinja2 are supported."
        )
