FROM nfcore/base:1.9
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN git clone https://github.com/AllanSSX/ElConcatenero.git /opt/ElConcatenero; \
    cd /opt/ElConcatenero; \
    chmod 777 *.py 
ENV PATH=/opt/conda/envs/pallor/bin:/opt/ElConcatenero:$PATH
