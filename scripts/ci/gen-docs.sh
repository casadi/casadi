#!/usr/bin/env bash
set -ex
cd "$(dirname "$0")"

output_folder="${1:-/tmp}"
mkdir -p "$output_folder"

# Function that builds the doxygen documentation and generates the
# coverage information.
# usage:    run_doxygen_coverage <branch-name> <output-directory>
function run_doxygen_coverage {
    repodir=$PWD/../..
    pushd ../../doxygen
    if [ "$1" = "main" ]; then dir="Doxygen"; else dir="$1/Doxygen"; fi
    # Remove the old documentation
    rm -rf "$2/$dir"
    mkdir -p "$2/$dir"
    # Tweak some Doxyfile verion numbers and output paths
    ( 
        cat Doxyfile;
        echo;
        echo "PROJECT_NUMBER = \"$1\"";
        echo "OUTPUT_DIRECTORY = \"$2\"";
        echo "HTML_OUTPUT = \"$dir\"";
        echo "GENERATE_LATEX = NO" 
    ) > tmp-Doxyfile
    # Generate the documentation
    doxygen tmp-Doxyfile
    rm -rf tmp-Doxyfile
    popd

    if [ "$1" = "main" ]; then dir="Coverage"; else dir="$1/Coverage"; fi
    pushd ../..
    rm -rf docs/Coverage build
    cmake -S. -Bbuild \
        -DALPAQA_WITH_COVERAGE=On \
        -DALPAQA_WITH_EXAMPLES=Off \
        -DALPAQA_WITH_CASADI=Off
    cmake --build build -j -t coverage
    mkdir -p "$2/$dir"
    rm -rf "$2/$dir"
    mv docs/Coverage "$2/$dir"
    popd
}

# Generate the documentation for the current branch
curr_branch=$(git branch --show-current)
if [ -n "$curr_branch" ]; then
    run_doxygen_coverage "$curr_branch" "$output_folder"
fi
# Generate the documentation for the current tag
if curr_tag=$(git describe --tags --exact-match); then
    run_doxygen_coverage "$curr_tag" "$output_folder"
fi

echo "Done generating documentation"

# Get all tags and branches for generating the index with links to docs for
# specific branches and versions:
git fetch
git fetch --tags

README="$output_folder/README.md"
echo "Documentation for" \
     "[**$GITHUB_REPOSITORY**](https://github.com/$GITHUB_REPOSITORY)." \
> "$README"
# Always have a link to main, it's at the root of the docs folder
echo -e '\n### Main Branch\n' >> "$README"
echo "- **main**  " >> "$README"
echo "  [Doxygen](Doxygen/index.html)"\
        "─ [Coverage](Coverage/index.html)" >> "$README"
# Find all tags with documentation (version numbers)
echo -e '\n### Tags and Releases\n' >> "$README"
git tag -l --sort=-creatordate \
 | while read tag; do
    index="$output_folder/$tag/Doxygen/index.html"
    if [ -e "$index" ]; then
        echo "- **$tag**  " >> "$README"
        echo "  [Doxygen]($tag/Doxygen/index.html)"\
             "─ [Coverage]($tag/Coverage/index.html)" >> "$README"
    else
        echo "tag $tag has no documentation"
    fi
done
# Find other branches (not version numbers)
echo -e '\n### Other Branches\n' >> "$README"
git branch -r --sort=-committerdate | cut -d / -f 2 \
 | while read branch
do
    index="$output_folder/$branch/Doxygen/index.html"
    if [ -e "$index" ]; then
        echo "- **$branch**  " >> "$README"
        echo "  [Doxygen]($branch/Doxygen/index.html)"\
             "─ [Coverage]($branch/Coverage/index.html)" >> "$README"
    else
        echo "branch $branch has no documentation"
    fi
done

echo -e "\n***\n" >> "$README"
# echo "<sup>Updated on $(date)</sup>" >> "$README"
cat > "$output_folder/_config.yml" << EOF
include:
  - "_modules"
  - "_sources"
  - "_static"
  - "_images"
EOF
