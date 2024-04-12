import tempfile
from pathlib import Path
from unittest.mock import patch

from snakemake.persistence import Persistence


class TestCleanupContainers:
    @patch("snakemake.deployment.singularity.Image", autospec=True, create=True)
    @patch("snakemake.dag.DAG", autospec=True, create=True)
    def test_one_unrequired_container_gets_removed(self, mock_dag, mock_img):
        snakecode = """rule all:
    input:
        "foo.txt"

rule foo:
    output:
        "foo.txt"
    container:
        "docker://quay.io/mbhall88/rasusa:0.7.0"
    shell:
        "rasusa --help &> {output}"
"""
        with tempfile.TemporaryDirectory() as tmpdirname:
            tmpdirpath = Path(tmpdirname)
            singularity_dir = tmpdirpath / ".snakemake/singularity"
            singularity_dir.mkdir(parents=True)
            snakefile = tmpdirpath / "Snakefile"
            snakefile.write_text(snakecode)

            unrequired_img_path = (
                singularity_dir / "32077ccc4d05977ef8b94ee6f74073fd.simg"
            )
            unrequired_img_path.touch()
            assert unrequired_img_path.exists()

            required_img_path = (
                singularity_dir / "0581af3d3099c1cc8cf0088c8efe1439.simg"
            )
            required_img_path.touch()
            assert required_img_path.exists()
            mock_img.path = str(required_img_path)
            container_imgs = {"docker://quay.io/mbhall88/rasusa:0.7.0": mock_img}
            mock_dag.container_imgs = container_imgs
            persistence = Persistence(dag=mock_dag)
            persistence.container_img_path = str(singularity_dir)

            persistence.cleanup_containers()

            assert required_img_path.exists()
            assert not unrequired_img_path.exists()
