# thesis
Repository for storing everything connected to my PhD thesis


## Git tokens
usethis::create_github_token()
gitcreds::gitcreds_set()

## Removing files you don't want Git-ted
git rm -r --cached .Rproj.user/**

## Downloading an entire folder from GitHub
https://download-directory.github.io

## Delete an un-pushed commit
git reset --soft HEAD~1 (keep changes in the working directory)
git reset --hard HEAD~1 (completely remove changes)

## The source function
https://www.earthdatascience.org/courses/earth-analytics/multispectral-remote-sensing-data/source-function-in-R