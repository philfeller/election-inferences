# Start with the election-inferences image

FROM ghcr.io/philfeller/election-inferences

# Install packages needed by the GitHub runner

RUN apt-get update && apt-get install -y curl \
    apt-transport-https \
    git \
    git-lfs \
    jq && \
    apt-get clean && apt-get autoremove

ARG GH_RUNNER_VERSION="2.331.0"

WORKDIR /actions-runner
COPY entrypoint.sh .

# Download the runner tarball, and extract the runner and its dependencies

# NOTE
# A GitHub token with full repo access is needed when a container is created and the
# runner is registered; it's passed as a secret but must be persisted to the image.
# Do not run the container in an environment where others would have access to it.
RUN --mount=type=secret,id=GH_TOKEN \
    curl -o actions.tar.gz --location "https://github.com/actions/runner/releases/download/v${GH_RUNNER_VERSION}/actions-runner-linux-x64-${GH_RUNNER_VERSION}.tar.gz" && \
    tar -zxf actions.tar.gz && \
    rm -f actions.tar.gz && \
    ./bin/installdependencies.sh && \
    chmod +x entrypoint.sh && \
    cp  /run/secrets/GH_TOKEN .gh_token

ENTRYPOINT ["/actions-runner/entrypoint.sh"]
