{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/23.05";
  };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };
    in
    {

      devShells.${system} = {
        pyshell = pkgs.mkShell {
          buildInputs = with pkgs; [
            python3
          ];
        };
        rshell = pkgs.mkShell {
          buildInputs = with pkgs; [
            R
          ];
        };
      };

    };
}
