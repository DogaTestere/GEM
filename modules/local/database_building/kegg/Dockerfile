FROM continuumio/miniconda3:latest

COPY modules/local/kegg/environment.yaml /tmp/environment.yaml

RUN conda env create -f /tmp/environment.yaml -n kegg_requests && \
    conda clean -afy

ENV PATH=/opt/conda/envs/kegg_requests/bin:$PATH

WORKDIR /work
CMD ["/bin/bash"]
