from snakemake.exceptions import WorkflowError
from snakemake.template_rendering import TemplateRenderer


class Jinja2Renderer(TemplateRenderer):
    def render(self):
        import jinja2

        try:
            with open(self.input_file, "r") as infile:
                template = jinja2.Template(infile.read())
            with open(self.output_file, "w") as outfile:
                outfile.write(template.render(**self.variables))
        except Exception as e:
            raise WorkflowError("Failed to render jinja2 template.", e)
