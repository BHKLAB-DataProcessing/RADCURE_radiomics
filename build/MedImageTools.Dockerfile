FROM python:3.8-slim-buster

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        && rm -rf /var/lib/apt/lists/*

RUN pip install med-imagetools==1.2.0.2 \
    && rm -rf /root/.cache/pip/*
