#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
repodir="$PWD"/../..
set -x
cd "$repodir"

mainbranch="main"
output_folder="${1:-/tmp}"
mkdir -p "$output_folder"

# Function that builds the doxygen documentation and generates the
# coverage information.
# usage:    run_doxygen_coverage <branch-name> <output-directory>
function run_doxygen_coverage {
    branch="$1"
    outdir="$2/$branch"
    htmldir="Doxygen"
    covdir="$2/$branch/Coverage"
    sphinxdir="$2/$branch/Sphinx"
    tmpdir="$repodir/tmp"
    # Prepare temporary folders
    mkdir -p "$tmpdir"
    echo '*' > "$tmpdir/.gitignore"
    # Remove the old documentation
    mkdir -p "$repodir/docs" "$outdir/$htmldir" "$covdir" "$sphinxdir"
    rm -rf "$outdir/$htmldir" "$covdir" "$sphinxdir"
    rm -rf "$repodir/docs/Coverage"

    # Tweak some Doxyfile verion numbers and output paths
    cat <<- EOF > "$tmpdir/tmp-Doxyfile"
	@INCLUDE = "$repodir/doxygen/Doxyfile"
	PROJECT_NUMBER = "$branch"
	OUTPUT_DIRECTORY = "$outdir"
	HTML_OUTPUT = "$htmldir"
	GENERATE_LATEX = NO
	EOF

    # Configure the project
    cmake -S. -B"$tmpdir/build" \
        -G "Ninja" \
        -DALPAQA_WITH_COVERAGE=On \
        -DALPAQA_WITH_TESTS=On \
        -DALPAQA_WITH_QUAD_PRECISION=On \
        -DALPAQA_WITH_PYTHON=Off \
        -DALPAQA_WITH_EXAMPLES=Off \
        -DALPAQA_WITH_CASADI=On \
        -DALPAQA_DOXYFILE="$tmpdir/tmp-Doxyfile"

    # Generate the Doxygen C++ documentation
    cmake --build "$tmpdir/build" -t docs

    # Tweak the settings for Doxygen for breathe
    cat <<- EOF > "$tmpdir/tmp-Doxyfile"
	@INCLUDE = "$repodir/doxygen/Doxyfile.breathe"
	PROJECT_NUMBER = "$branch"
	OUTPUT_DIRECTORY = "$tmpdir"
	XML_OUTPUT = "xml"
	GENERATE_LATEX = NO
	EOF

    # Generate the Sphinx Python & C++ documentation
    cmake --build "$tmpdir/build" -t docs
    sphinx-build -b doctest -j auto -D "breathe_projects.alpaqa=$tmpdir/xml" \
        doxygen/sphinx/source "$sphinxdir" \
    || echo -e "\n##################\n# DOCTEST FAILED #\n##################\n"
    sphinx-build -b html -j auto -D "breathe_projects.alpaqa=$tmpdir/xml" \
        doxygen/sphinx/source "$sphinxdir"

    # Generate coverage report
    cmake --build "$tmpdir/build" -j -t coverage
    mv docs/Coverage "$covdir"

    # Cleanup
    rm -f tmp-Doxyfile
}

# Generate the documentation for the current branch
git fetch ||:
curr_branch=$(git branch --show-current)
if [ -n "$curr_branch" ]; then
    run_doxygen_coverage "$curr_branch" "$output_folder"
elif [ -n "$CI_COMMIT_BRANCH" ]; then
    run_doxygen_coverage "$CI_COMMIT_BRANCH" "$output_folder"
fi
# Generate the documentation for the current tag
git fetch --tags ||:
if curr_tag=$(git describe --tags --exact-match); then
    run_doxygen_coverage "$curr_tag" "$output_folder"
fi

set +x

echo "Done generating documentation"
