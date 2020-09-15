__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2020, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os
import re
import requests


class DockerContainer:
    """A Docker Container is a thin wrapper to a Docker Registry API to
    primarily determine if a container exists to be used for a workflow
    """

    def __init__(self, name):
        self.name = name
        self.uri = None
        if self.exists():
            self.uri = "docker://%s" % self.name

    @property
    def version(self):
        """Derive a container verison. Will by default be latest (the tag)
        but can be a different tag or full hash.
        """
        # Does the container have a version (tag or hash?)
        version = "latest"
        if "@" in self.name:
            container_name, version = self.name.split("@", 1)
        elif ":" in self.name:
            container_name, version = self.name.split(":", 1)
        return version

    @property
    def registry(self):
        """Return the container registry. We default to Docker Hub unless
        another (quay.io) is defined
        """
        registry = "https://index.docker.io/v2"
        if "quay.io" in self.name:
            registry = "https://quay.io/api/v1"
        return registry

    @property
    def repository(self):
        """The repository is the username/reponame combination, without
        the registry. E.g., quay.io/name/repo --> name/repo
        """
        container_name = self.name.replace("quay.io/", "", 1)
        return re.split("(:|@)", container_name, 1)[0]

    def exists(self):
        """Return True if the container exists under CONTAINER_PREFIX
        (quay.io/snakemake-wrappers) otherwise return False
        """
        # Cut out early if given no Name, possibly creating instance for later use.
        if not self.name:
            return

        # If the user provides a docker uri, ensure we remove it first
        self.name = self.name.replace("docker://", "", 1)
        if "quay.io" in self.registry:
            return self.exists_quay()
        return self.exists_docker()

    def exists_quay(self):
        """Quay will return a 200 response if a repository, tag combination exists,
        and return a json object with "images" -> images. We only care if the container
        exists, which we can determine based on the response.status_code.
        """
        # get /api/v1/repository/{repository}/tag/{tag}/images
        url = "%s/repository/%s/tag/%s/images" % (
            self.registry,
            self.repository,
            self.version,
        )

        response = requests.get(url)
        return response.status_code == 200

    def exists_docker(self):
        """Determine if an image exists for docker hub. We are able to query for
        an image manifest (referenced by tag) directly after getting an auth
        token.
        """
        # Request for the manifest
        # GET /v2/<name>/manifests/<reference>
        base = "%s/%s/manifests/%s" % (self.registry, self.repository, self.version)
        headers = {"Docker-Distribution-API-Version": "registry/2.0"}
        response = requests.get(base, headers=headers)

        # Repository does not exist
        if response.status_code == 404:
            return False

        # Authenticate
        elif response.status_code == 401 and "Www-Authenticate" in response.headers:
            challenge = response.headers["Www-Authenticate"]
            token_url = self.get_token_url(challenge, expires_in=900)

        # Get token to authenticate
        data = requests.get(token_url, headers=headers).json()
        if "access_token" in data:
            access_token = data["access_token"]
        else:
            access_token = data["token"]
        token = {"Authorization": "Bearer %s" % access_token}
        headers.update(token)
        response = requests.get(base, headers=headers)
        return response.status_code == 200

    def get_token_url(self, challenge, expires_in, sort_query_params=False):
        """Build token URL from authentication challenge"""
        params = parse_bearer_challenge(challenge)

        if params and "realm" in params:
            realm = params.pop("realm")
            params["expires_in"] = expires_in

        items = params.items()
        if sort_query_params:
            items = sorted(params.items())

        query_fragment = "&".join(["%s=%s" % (k, v) for k, v in items])
        return "{realm}?{query_fragment}".format(
            realm=realm, query_fragment=query_fragment
        )


def parse_bearer_challenge(challenge):
    """Parse a bearer challenge from a 'Www-Authenticate' header into a
    dictionary of params. Supports realm, service, and scope bearer params
    appearing in any order and assumes param values are double quoted.
    """
    param_re = '((?:realm|service|scope)=".+?")'
    regexp = "^Bearer\s+{param},?{param}?,?{param}?".format(param=param_re)
    match = re.match(regexp, challenge)
    params = {}
    if match:
        for group in match.groups():
            if group:
                param, value = group.split("=", 1)
                params[param] = value.strip('"')
    return params
