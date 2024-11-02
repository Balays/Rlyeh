#!/bin/bash

# Usage message
usage() {
  echo "Usage: $0 bamfile tag"
  echo "  bamfile : Path to the BAM file."
  echo "  tag     : The tag to count values for (e.g., 'fn')."
  exit 1
}

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
  usage
fi

BAMFILE=$1
TAG=$2

# Check if the BAM file exists
if [ ! -f "$BAMFILE" ]; then
  echo "Error: BAM file $BAMFILE not found."
  exit 1
fi

# Extract the specified tag values and count occurrences
samtools view "$BAMFILE" | \
awk -v tag="$TAG" '{
  for (i=12; i<=NF; i++) {
    split($i, a, ":")
    if (a[1] == tag) {
      values[a[3]]++
    }
  }
} END {
  for (value in values) {
    print value, values[value]
  }
}' | sort -k2 -nr
