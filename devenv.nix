{ pkgs, lib, config, inputs, ... }:

let
  # --- Pinned Nextflow -------------------------------------------------
  nextflowVersion = "26.04.6";

  nextflowPinned = pkgs.stdenv.mkDerivation {
    pname = "nextflow";
    version = nextflowVersion;

    src = pkgs.fetchurl {
      url = "https://github.com/nextflow-io/nextflow/releases/download/v${nextflowVersion}/nextflow-${nextflowVersion}-dist";
      sha256 = "sha256-GCpjx0B04tx5Vv+jyM1Z3pUu0sRDlOIfr14XNrlFREw=";
    };

    dontUnpack = true;
    nativeBuildInputs = [ pkgs.makeWrapper ];

    installPhase = ''
      mkdir -p $out/bin
      install -m755 $src $out/bin/nextflow
      wrapProgram $out/bin/nextflow --prefix PATH : ${lib.makeBinPath [ pkgs.jdk21 ]}
    '';
  };
in

{
  env.GREET = "Metabolic Project";

  packages = [
    nextflowPinned
    pkgs.docker-client 
    pkgs.micromamba      
    pkgs.python312       
  ];

  scripts.hello.exec = ''
    echo hello from $GREET
  '';

  enterShell = ''
    if [ ! -d .venv ]; then
      python3 -m venv .venv
    fi
    source .venv/bin/activate

    if [ -f requirements.txt ]; then
      pip install -q -r requirements.txt
    fi

    echo hello
    echo "  nextflow: $(nextflow -v 2>/dev/null | head -n1)"
    echo "  nf-core:  $(nf-core --version 2>/dev/null)"
    echo "  conda:    $(micromamba --version 2>/dev/null)"
    echo "  docker:   $(docker --version 2>/dev/null || echo 'daemon not reachable — see README')"
  '';
}
