from snakemake.exceptions import WorkflowError
from snakemake.template_rendering import TemplateRenderer


class YteRenderer(TemplateRenderer):
    def render(self):
        import yte

        try:
            with (
                open(self.output_file, "w") as outfile,
                open(self.input_file, "r") as infile,
            ):
                yte.process_yaml(infile, outfile=outfile, variables=self.variables)
        except Exception as e:
            raise WorkflowError("Failed to render yte template.", e)
