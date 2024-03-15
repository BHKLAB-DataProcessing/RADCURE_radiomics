# channels:
#   - defaults
#   - anaconda
#   - conda-forge
#   - bioconda
# dependencies:
#   - numpy=1.26.4
#   - python=3.9.18
#   - pip=23.3.1
#   - pip:
#       - readii==1.1.3

FROM python:3.9-slim-buster

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        && rm -rf /var/lib/apt/lists/*

# change python version
RUN pip install readii==1.1.3
RUN pip install numpy==1.26.4
