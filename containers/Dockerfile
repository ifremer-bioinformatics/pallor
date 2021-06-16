FROM nfcore/base:1.9
COPY environment.yml /
RUN conda env create -f /environment.yml; \
    conda clean -a
WORKDIR /opt
RUN git clone https://github.com/ODiogoSilva/ElConcatenero.git; \
    chmod +x ElConcatenero/*.py
ENV PATH=/opt/conda/envs/pallor/bin:/opt/ElConcatenero:$PATH
