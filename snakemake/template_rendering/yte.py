from snakemake.exceptions import WorkflowError
from snakemake.template_rendering import TemplateRenderer


class YteRenderer(TemplateRenderer):
    def render(self):
        import yte

        try:
            yte.process_yaml(
                self.input_file, outfile=self.output_file, variables=self.variables
            )
        except Exception as e:
            raise WorkflowError("Failed to render yte template.", e)
