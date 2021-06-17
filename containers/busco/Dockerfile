FROM nfcore/base:1.9
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH=/opt/conda/envs/busco/bin:$PATH
