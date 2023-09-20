

class NixFlakeEnv:

    """Nix Flake Env"""

    def __init__(self, nix_flake):
        self.nix_flake = nix_flake

    def shellcmd(self, cmd):
        return f"nix develop {self.nix_flake} --command bash -c '{cmd}'"
