#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
repodir="$PWD"/../..
set -x
cd "$repodir"

mainbranch="develop"
output_folder="${1:-/tmp}"
mkdir -p "$output_folder"

# Function that builds the doxygen documentation and generates the
# coverage information.
# usage:    run_doxygen_coverage <branch-name> <output-directory>
function run_doxygen_coverage {
    branch="$1"
    if [ "$branch" = "$mainbranch" ]; then
        outdir="$2"
    else
        outdir="$2/$branch"
    fi
    htmldir="Doxygen"
    covdir="$outdir/Coverage"
    sphinxdir="$outdir/Sphinx"
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

    # Install dependencies
    conan build . -pr scripts/ci/profiles/x86_64-bionic-linux-gnu.profile \
        --build=missing \
        -s \&:build_type=Debug \
        -c \*:tools.build:skip_test=True \
        -c \&:tools.build:skip_test=False \
        -c \&:tools.build.cross_building:can_run=true \
        -c \&:tools.cmake.cmaketoolchain:generator=Ninja \
        -c \&:tools.cmake.cmaketoolchain:extra_variables="{\"ALPAQA_DOXYFILE\": \"$tmpdir/tmp-Doxyfile\"}" \
        -c \&:tools.cmake.cmake_layout:build_folder_vars="['const.docs']" \
        -o \&:with_coverage=True \
        -o \&:with_quad_precision=False \
        -o \&:with_examples=False \
        -o \&:with_casadi=True \
        -o \&:with_cutest=True

    # Generate the Doxygen C++ documentation
    cmake --build --preset conan-docs-debug -t docs

    # Need to support old links to modules.html
    [ -e "$outdir/$htmldir/modules.html" ] || \
        cp "$outdir/$htmldir/topics.html" "$outdir/$htmldir/modules.html"

    # Tweak the settings for Doxygen for breathe
    cat <<- EOF > "$tmpdir/tmp-Doxyfile"
	@INCLUDE = "$repodir/doxygen/Doxyfile.breathe"
	PROJECT_NUMBER = "$branch"
	OUTPUT_DIRECTORY = "$tmpdir"
	XML_OUTPUT = "xml"
	GENERATE_LATEX = NO
	EOF

    # Generate the Sphinx Python & C++ documentation
    cmake --build --preset conan-docs-debug -t docs
    sphinx-build -b doctest -j auto -D "breathe_projects.alpaqa=$tmpdir/xml" \
        doxygen/sphinx/source "$sphinxdir" \
    || echo -e "\n##################\n# DOCTEST FAILED #\n##################\n"
    sphinx-build -b html -j auto -D "breathe_projects.alpaqa=$tmpdir/xml" \
        doxygen/sphinx/source "$sphinxdir"

    # Generate coverage report
    cmake --build --preset conan-docs-debug -t coverage
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
