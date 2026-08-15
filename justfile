#!/usr/bin/env just --justfile

latest-tag := "0.2.0"
# Matches the crate and binary name. Was `dlejeune/pipeline-utils-rs`, which matched
# nothing else in the repo.
image-name := "dlejeune/rusty-metal"

default:
    just --list

# Builds a docker image with the most recent git tag
[group('docker')]
build-docker:
    sudo docker build -t {{ image-name }}:{{ latest-tag }} .

# Pushes the docker image with the most recent git tag to dockerhub
[group('docker')]
push-docker:
    sudo docker push {{ image-name }}:{{ latest-tag }}

# Runs an interactive docker container with the current wd mounted at /data
[group('docker')]
run-docker-it tag=latest-tag:
    sudo docker run --rm -it -v ./:/data {{ image-name }}:{{ tag }} bash

# Builds and pushed the most recently tagged branch in a docker container
[group('docker')]
docker: build-docker push-docker

build:
    cargo build

run *args="":
    cargo run -- {{ args }}

# Runs the test suite
test *args="":
    cargo test {{ args }}

# Formatting, lints and tests — what CI should run. All three are clean as of 0.2.0.
check:
    cargo fmt --check
    cargo clippy -- -D warnings
    cargo test

# Formats the tree in place. Run this before `check` if it complains about formatting.
fmt:
    cargo fmt
