FROM rust:1 as build-env
WORKDIR /app
COPY . /app
RUN rustup target add x86_64-unknown-linux-musl
RUN apt update && \
    apt install -y musl-tools musl-dev && \
    update-ca-certificates
RUN cargo build --target x86_64-unknown-linux-musl --release

FROM alpine:latest as release
COPY --from=build-env /app/target/x86_64-unknown-linux-musl/release/seq2c-rs /bin/seq2c-rs
CMD ["./bin/seq2c-rs"]
