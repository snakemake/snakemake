#!/bin/bash

# Extract version from PR tag passed as environment variable
if [ -z "${PR_TITLE}" ]; then
    echo "Error: 'PR_TITLE' environment variable is not set."
    exit 1
fi

# Check if this is a release PR
if [[ ! "${PR_TITLE}" =~ ^chore\(main\):\ release ]]; then
    echo "Not a release PR, skipping Mastodon post"
    exit 0
fi

# Extract version (everything after "release ")
if [[ "${PR_TITLE}" =~ [Rr]elease[[:space:]]+([0-9]+\.[0-9]+\.[0-9]+) ]]; then
    version="${BASH_REMATCH[1]}"
else
    echo "Error: Could not extract version number from PR title"
    exit 1
fi

# Validate version format
if ! [[ $version =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "Error: Invalid version format in PR title: $version"
    exit 1
fi

# Construct changelog URL with proper quoting
changelog="https://github.com/snakemake/snakemake-executor-plugin-slurm/releases/tag/v${version}"

# Maximum character limit for Mastodon posts (on Fediscience: 1500 characters)
MAX_TOOT_LENGTH=1500

read -d '\n' message << EndOfText
BEEP, BEEP - I  am your friendly #Snakemake release announcement bot.

There is a new release of Snakemake itself. Its version is now '${version}'!

See ${changelog//\'/\\\'} for details.

Give us some time and you will automatically find the plugin on #Bioconda and #Pypi.

If you want to discuss the release you will find the maintainer here on Mastodon! - @johanneskoester@fosstodon.org

If you find any issues, please report them on https://github.com/snakemake/snakemake/issues .

#Snakemake #ReproducibleComputing #ReproducibleResearch #OpenScience
EndOfText

# Validate message length
if [ ${#message} -gt $MAX_TOOT_LENGTH ]; then
    echo "Error: Message exceeds Fediscience's character limit"
    exit 1
fi

# Validate Mastodon token
if [ -z "${MASTODONBOT}" ]; then
    echo "Error: MASTODONBOT secret is not set"
    exit 1
fi

# Send post to Mastodon with proper quoting and error handling
response=$(curl -s -w "\n%{http_code}" -X POST \
    -H "Authorization: Bearer ${MASTODONBOT}" \
    -F "status=${message}" \
    "https://fediscience.org/api/v1/statuses")

status_code=$(echo "$response" | tail -n1)
response_body=$(echo "$response" | sed '$d')

if [ "$status_code" -ne 200 ]; then
    echo "Error: Failed to post to Mastodon (HTTP ${status_code})"
    echo "Response: ${response_body}"
    exit 1
fi

echo "Successfully posted to Mastodon"
