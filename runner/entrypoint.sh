#!/usr/bin/env bash
set -eEuo pipefail

ACTIONS_RUNNER_INPUT_NAME=$HOSTNAME
export RUNNER_ALLOW_RUNASROOT=1

# Generate a token for registering the runner; authorization requires a token with full repo access
# This token is stored in action secrets and passed to the build as a secret, but it must be
# persisted to the image to be available at run time.
TOKEN="$(curl -sS --request POST --url https://api.github.com/repos/philfeller/election-inferences/actions/runners/registration-token --header "authorization: Bearer $( cat .gh_token)" --header 'content-type: application/json' | jq -r .token)"

# Start the runner, which registers itself as self-hosted
/actions-runner/config.sh --url https://github.com/philfeller/election-inferences --token "$TOKEN" --unattended --labels election-inferences
/actions-runner/bin/runsvc.sh
