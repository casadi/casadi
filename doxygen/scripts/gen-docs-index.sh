#!/usr/bin/env bash
set -e

output_folder="$(realpath "${1:-/tmp}")"
output_type="${2:-md}"

cd "$(dirname "$0")"
repodir="$PWD"/../..
set -x
cd "$repodir"

mainbranch="main"

set +x

# Get all tags and branches for generating the index with links to docs for
# specific branches and versions:
git fetch
git fetch --tags

function write_readme_item_item {
    type="$1"
    name="$2"
    path="$3"
    itemname="${4:-$3}"
    if [ -e "$output_folder/$name/$path/index.html" ]; then
        echo "  [$itemname]($name/$path/)  \n"
    else
        echo "$type $name has no documentation for $itemname" >&2
    fi
}

function write_readme_item {
    type="$1"
    name="${2//\//-}"
    README="$3"
    local result=""
    result+="$(write_readme_item_item "$type" "$name" Doxygen)"
    result+="$(write_readme_item_item "$type" "$name" Sphinx)"
    # Write result to file
    if [ -n "$result" ]; then
        echo "- **$name**  " >> "$README"
        echo -e "$result" >> "$README"
    fi
}

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
    write_readme_item "main branch" "$mainbranch" "$README"


    # Find all tags with documentation (version numbers)
    echo -e '\n### Tags and releases\n' >> "$README"
    git tag -l --sort=-creatordate \
    | while read tag
    do
        write_readme_item "tag" "$tag" "$README"
    done

    # Find other branches (not version numbers)
    echo -e '\n### Other branches\n' >> "$README"
    git branch -r --sort=-committerdate | cut -d / -f 2 \
    | while read branch
    do
        if [ "$branch" = "$mainbranch" ]; then
            : # skip the main branch
        else
            write_readme_item "branch" "$branch" "$README"
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

function write_index_item_item {
    type="$1"
    name="$2"
    path="$3"
    itemname="${4:-$3}"
    if [ -e "$output_folder/$name/$path/index.html" ]; then
        echo "<br>\n  <a href=\"$name/$path/\">$itemname</a>"
    else
        echo "$type $name has no documentation for $itemname" >&2
    fi
}

function write_index_item {
    type="$1"
    name="${2//\//-}"
    README="$3"
    local result=""
    result+="$(write_index_item_item "$type" "$name" Doxygen)"
    result+="$(write_index_item_item "$type" "$name" Sphinx)"
    # Write result to file
    if [ -n "$result" ]; then
        echo -n "<li><b>$name</b>" >> "$README"
        echo -e "$result" >> "$README"
        echo -e "\n</li>" >> "$README"
    fi
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
    write_index_item "main branch" "$mainbranch" "$README"
    echo -e '</ul>' >> "$README"

    # Find all tags with documentation (version numbers)
    echo -e '\n<h3>Tags and releases</h3>\n' >> "$README"
    echo -e '<ul>' >> "$README"
    git tag -l --sort=-creatordate \
    | while read tag
    do
        write_index_item "tag" "$tag" "$README"
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
        else
            write_index_item "branch" "$branch" "$README"
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