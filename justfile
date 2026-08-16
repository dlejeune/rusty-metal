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

# Builds the API documentation, with warnings fatal.
#
# `cargo rustdoc` rather than `cargo doc` because it is the form that takes rustdoc
# flags directly, so `-D warnings` needs no RUSTDOCFLAGS in the environment. Private
# items are documented without asking: this is a binary target, and every module in it
# is private to `main.rs`, so the default for libraries would produce an empty page.
doc *args="":
    cargo rustdoc --bin rusty-metal -- -D warnings {{ args }}

# Builds the documentation and opens it in a browser.
doc-open:
    cargo doc --no-deps --open

# Formatting, lints, docs and tests — what CI runs. All four are clean as of 0.2.0.
#
# `doc` is in here because a doc comment can rot in ways nothing else catches: a
# `[link]` to an item that was renamed, or made `#[cfg(test)]`, still compiles and still
# passes the tests. That is exactly how the reference to `Element` in `homology.rs`
# survived Stage 8.
check:
    cargo fmt --check
    cargo clippy -- -D warnings
    just doc
    cargo test

# Formats the tree in place. Run this before `check` if it complains about formatting.
fmt:
    cargo fmt
