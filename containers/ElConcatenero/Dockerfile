FROM quay.io/biocontainers/biopython:1.78
COPY ElConcatenero /opt/ElConcatenero
RUN cd /opt/ElConcatenero; \
    chmod 777 *.py
ENV PATH=/opt/ElConcatenero:$PATH
