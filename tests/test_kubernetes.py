import os
import sys
import uuid

sys.path.insert(0, os.path.dirname(__file__))

from common import *


@pytest.fixture(scope="module")
def kubernetes_cluster():
    class Cluster:
        def __init__(self):
            self.cluster = "t-{}".format(uuid.uuid4())
            self.bucket_name = self.cluster

            shell(
                """
                gcloud container clusters create {self.cluster} --num-nodes 3 --scopes storage-rw --zone us-central1-a --machine-type n1-standard-2
                gcloud container clusters get-credentials {self.cluster} --zone us-central1-a
                gsutil mb gs://{self.bucket_name}
                """
            )

        def delete(self):
            shell(
                """
                gcloud container clusters delete {self.cluster} --zone us-central1-a --quiet || true
                gsutil rm -r gs://{self.bucket_name} || true
                """
            )

        def run(self, test="test_kubernetes", **kwargs):
            try:
                run(
                    dpath(test),
                    kubernetes="default",
                    default_remote_provider="GS",
                    default_remote_prefix=self.bucket_name,
                    no_tmpdir=True,
                    **kwargs
                )
            except Exception as e:
                shell(
                    "for p in `kubectl get pods | grep ^snakejob- | cut -f 1 -d ' '`; do kubectl logs $p; done"
                )
                raise e

        def reset(self):
            print("Resetting bucket...", file=sys.stderr)
            shell("gsutil -m rm -r gs://{self.bucket_name}/* || true")

    cluster = Cluster()
    yield cluster
    cluster.delete()


@gcloud
def test_kubernetes_plain(kubernetes_cluster):
    kubernetes_cluster.reset()
    kubernetes_cluster.run()


@gcloud
@pytest.mark.skip(reason="need a faster cloud compute instance to run this")
def test_kubernetes_conda(kubernetes_cluster):
    kubernetes_cluster.reset()
    kubernetes_cluster.run(use_conda=True)


@gcloud
@pytest.mark.skip(reason="need a faster cloud compute instance to run this")
def test_kubernetes_singularity(kubernetes_cluster):
    kubernetes_cluster.reset()
    kubernetes_cluster.run(use_singularity=True)


@gcloud
@pytest.mark.skip(reason="need a faster cloud compute instance to run this")
def test_kubernetes_conda_singularity(kubernetes_cluster):
    kubernetes_cluster.reset()
    kubernetes_cluster.run(use_singularity=True, use_conda=True)


@gcloud()
@pytest.mark.skip(reason="need a faster cloud compute instance to run this")
def test_issue1041(kubernetes_cluster):
    kubernetes_cluster.reset()
    kubernetes_cluster.run(test="test_issue1041")
