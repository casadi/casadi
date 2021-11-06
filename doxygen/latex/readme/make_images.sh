#!/usr/bin/env bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

inkscape problem.pdf --pdf-poppler --export-text-to-path --export-plain-svg -o ../../images/problem.svg
