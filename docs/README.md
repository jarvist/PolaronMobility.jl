# Documentation

Note to self: Built with Documenter.jl (extracting doc-strings) then Python's mkdocs.

WILL NOW (May 2022) BUILD AND SELF DEPLOY VIA GITHUB ACTION

## Manual instructions

julia make.jl # pull in latest doc strings

mkdocs serve # local live-serve edits

mkdocs gh-deploy # Push to https://jarvist.github.io/PolaronMobility.jl/

