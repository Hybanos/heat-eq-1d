FROM ubuntu:20.04

RUN apt-get update -y && \
  apt-get install -y gcc libblas-dev liblapacke-dev make vim valgrind gdb

WORKDIR /app
# COPY . .

ENV HOSTNAME=docker
# RUN make all && make run