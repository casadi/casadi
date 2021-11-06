#!/usr/bin/env bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

inkscape alpaqa_favicon.pdf --pdf-poppler --export-text-to-path --export-plain-svg -o alpaqa_favicon.svg
inkscape alpaqa_logo.pdf --pdf-poppler --export-text-to-path --export-plain-svg -o alpaqa_logo.svg
inkscape alpaqa_favicon.pdf --pdf-poppler --export-type png -w 256 -h 256
