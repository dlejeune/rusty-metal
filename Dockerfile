# syntax=docker/dockerfile:1
#
# The image the justfile's `build-docker` / `push-docker` / `run-docker-it` recipes use.
#
# Two stages: a toolchain image that compiles, and a bare Debian that carries only the
# binary. The build stage is ~1.5 GB and none of it is needed at runtime.

# Pinned rather than `rust:1-slim`, so a rebuild six months from now produces the same
# binary. Edition 2024 needs 1.85 or newer; 1.94 is what the tree is developed against.
FROM rust:1.94-slim-bookworm AS build

WORKDIR /src
COPY . .

# BuildKit cache mounts for the registry and the target directory, so an edit to one
# source file does not recompile every dependency. The binary has to be copied out
# inside the same RUN: a cache mount is not part of the resulting layer, so anything
# left in /src/target vanishes when this instruction finishes.
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/src/target \
    cargo build --release --locked \
    && cp target/release/rusty-metal /usr/local/bin/rusty-metal

# `--locked` because Cargo.lock is committed: a build that silently picked up newer
# dependencies than the tested tree would produce an image nobody can reproduce.

FROM debian:bookworm-slim

# Same Debian release as the build stage. The binary links against the system glibc, so
# these two have to agree — a newer builder against an older runtime fails at startup
# with a GLIBC_2.xx not found that is tedious to diagnose.

COPY --from=build /usr/local/bin/rusty-metal /usr/local/bin/rusty-metal

# Where `just run-docker-it` mounts the host's working directory.
WORKDIR /data

# No ENTRYPOINT, on purpose. `just run-docker-it` runs `docker run ... <image> bash`, and
# an entrypoint of `rusty-metal` would feed it "bash" as a positional argument instead of
# starting a shell. The cost is that the binary has to be named explicitly:
#
#     docker run --rm -v ./:/data dlejeune/rusty-metal:0.2.0 \
#         rusty-metal -o /data/dist.csv /data/a.fasta /data/b.fasta
#
CMD ["rusty-metal", "--help"]

# Runs as root. The tool reads and writes files in a bind-mounted directory, and a
# baked-in non-root UID cannot write to a host directory owned by anyone else, which
# would leave every user passing `--user $(id -u):$(id -g)`. There is no network listener
# here and nothing persistent in the image.
