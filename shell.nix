{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = with pkgs; [
    root
    boost
    gnumake
    gcc
  ];

}