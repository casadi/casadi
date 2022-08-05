#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
repodir="$PWD"/../..
set -x
cd "$repodir"

mainbranch="main"
output_folder="${1:-/tmp}"
output_type="${2:-md}"
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
    cat <<- EOF > tmp-Doxyfile
	@INCLUDE = "$repodir/doxygen/Doxyfile"
	PROJECT_NUMBER = "$branch"
	OUTPUT_DIRECTORY = "$outdir"
	HTML_OUTPUT = "$htmldir"
	GENERATE_LATEX = NO
	EOF

    # Configure the project
    cmake -S. -B"$tmpdir/build" \
        -DALPAQA_WITH_COVERAGE=On \
        -DALPAQA_WITH_TESTS=On \
        -DALPAQA_WITH_QUAD_PRECISION=On \
        -DALPAQA_WITH_PYTHON=Off \
        -DALPAQA_WITH_EXAMPLES=Off \
        -DALPAQA_WITH_CASADI=Off \
        -DALPAQA_DOXYFILE="$repodir/tmp-Doxyfile"

    # Generate the Doxygen C++ documentation
    cmake --build "$tmpdir/build" -t docs

    # Generate the Sphinx Python & C++ documentation
    cmake "$tmpdir/build" -DALPAQA_DOXYFILE="$repodir/doxygen/Doxyfile.breathe"
    cmake --build "$tmpdir/build" -t docs
    sphinx-build -b html -j auto doxygen/sphinx/source "$sphinxdir"

    # Generate coverage report
    cmake --build "$tmpdir/build" -j -t coverage
    mv docs/Coverage "$covdir"

    # Cleanup
    rm -f tmp-Doxyfile
}

# Get all tags and branches for generating the index with links to docs for
# specific branches and versions:
git fetch
git fetch --tags

# Generate the documentation for the current branch
curr_branch=$(git branch --show-current)
if [ -n "$curr_branch" ]; then
    run_doxygen_coverage "$curr_branch" "$output_folder"
elif [ -n "$CI_COMMIT_BRANCH" ]; then
    run_doxygen_coverage "$CI_COMMIT_BRANCH" "$output_folder"
fi
# Generate the documentation for the current tag
if curr_tag=$(git describe --tags --exact-match); then
    run_doxygen_coverage "$curr_tag" "$output_folder"
fi

set +x

echo "Done generating documentation"

function write_readme {

    README="$output_folder/README.md"
    if [ -n "$GITHUB_REPOSITORY" ]; then
        echo "Documentation for" \
            "[**$GITHUB_REPOSITORY**](https://github.com/$GITHUB_REPOSITORY)." \
        > "$README"
    elif [ -n "$CI_PROJECT_URL" ]; then
        echo "Documentation for" \
            "[**$CI_PROJECT_PATH**]($CI_PROJECT_URL)." \
        > "$README"
    else
        echo "Documentation." \
        > "$README"
    fi


    # Always have a link to main, it's at the root of the docs folder
    echo -e '\n### Main branch\n' >> "$README"
    echo "- **$mainbranch**  " >> "$README"
    echo "  [Doxygen]($mainbranch/Doxygen/index.html)" >> "$README"

    # Find all tags with documentation (version numbers)
    echo -e '\n### Tags and releases\n' >> "$README"
    git tag -l --sort=-creatordate \
    | while read tag
    do
        index="$output_folder/$tag/Doxygen/index.html"
        if [ -e "$index" ]; then
            echo "- **$tag**  " >> "$README"
            echo "  [Doxygen]($tag/Doxygen/index.html)" >> "$README"
        else
            echo "tag $tag has no documentation"
        fi
    done

    # Find other branches (not version numbers)
    echo -e '\n### Other branches\n' >> "$README"
    git branch -r --sort=-committerdate | cut -d / -f 2 \
    | while read branch
    do
        index="$output_folder/$branch/Doxygen/index.html"
        if [ "$branch" = "$mainbranch" ]; then
            : # skip the main branch
        elif [ -e "$index" ]; then
            echo "- **$branch**  " >> "$README"
            echo "  [Doxygen]($branch/Doxygen/index.html)" >> "$README"
        else
            echo "branch $branch has no documentation"
        fi
    done

    echo -e "\n***\n" >> "$README"
    # echo "<sup>Updated on $(date)</sup>" >> "$README"
    cat > "$output_folder/_config.yml" <<- EOF
	include:
	- "_modules"
	- "_sources"
	- "_static"
	- "_images"
	EOF

}

function write_index {

    README="$output_folder/index.html"
    cat "$repodir/doxygen/custom/index-header.html" > "$README"
    if [ -n "$GITHUB_REPOSITORY" ]; then
        echo "<p>Documentation for" \
            "<a href=\"https://github.com/$GITHUB_REPOSITORY\"><b>$GITHUB_REPOSITORY</b></a>." \
            "</p>" \
        >> "$README"
    elif [ -n "$CI_PROJECT_URL" ]; then
        echo "<p>Documentation for" \
            "<a href=\"$CI_PROJECT_URL\"><b>$CI_PROJECT_PATH</b></a>.</p>" \
        >> "$README"
    else
        echo "<p>Documentation.</p>" \
        >> "$README"
    fi


    # Always have a link to main, it's at the root of the docs folder
    echo -e '\n<h3>Main branch</h3>\n' >> "$README"
    echo -e '<ul>' >> "$README"
    echo "<li><b>$mainbranch</b><br>" >> "$README"
    echo "<a href=\"$mainbranch/Doxygen/index.html\">Doxygen</a></li>" >> "$README"
    echo -e '</ul>' >> "$README"

    # Find all tags with documentation (version numbers)
    echo -e '\n<h3>Tags and releases</h3>\n' >> "$README"
    echo -e '<ul>' >> "$README"
    git tag -l --sort=-creatordate \
    | while read tag
    do
        index="$output_folder/$tag/Doxygen/index.html"
        if [ -e "$index" ]; then
            echo "<li><b>$tag</b><br>" >> "$README"
            echo "<a href=\"$tag/Doxygen/index.html\">Doxygen</a></li>" >> "$README"
        else
            echo "tag $tag has no documentation"
        fi
    done
    echo -e '</ul>' >> "$README"

    # Find other branches (not version numbers)
    echo -e '\n<h3>Other branches</h3>\n' >> "$README"
    echo -e '<ul>' >> "$README"
    git branch -r --sort=-committerdate | cut -d / -f 2 \
    | while read branch
    do
        index="$output_folder/$branch/Doxygen/index.html"
        if [ "$branch" = "$mainbranch" ]; then
            : # skip the main branch
        elif [ -e "$index" ]; then
            echo "<li><b>$branch</b><br>" >> "$README"
            echo "<a href=\"$branch/Doxygen/index.html\">Doxygen</a></li>" >> "$README"
        else
            echo "branch $branch has no documentation"
        fi
    done
    echo -e '</ul>' >> "$README"

    cat "$repodir/doxygen/custom/index-footer.html" >> "$README"

}

case "$output_type" in 
    md) write_readme ;;
    html) write_index ;;
    *) echo "Invalid output type $output_type"; exit 1 ;;
esac
