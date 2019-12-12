#!/usr/bin/env bash

# Get all occurrences of use m_... in source files
# deps="$(grep -r -e "use m_" --include \*.f90)"
deps="$(find . -maxdepth 1 -iname "m_*.f90" -type f -exec grep -H "use m_" {} \;)"

# Remove comments
deps=$(echo "$deps" | sed 's/!.*$//')

# Remove 'only: ...'
deps=$(echo "$deps" | sed 's/, *only.*$//')

# Remove m_af_...
deps=$(echo "$deps" | sed 's/use m_af_.*$//')

# Remove m_particle_core and related
deps=$(echo "$deps" | sed 's/use m_particle_core.*$//')
deps=$(echo "$deps" | sed 's/use m_gas.*$//')
deps=$(echo "$deps" | sed 's/use m_cross_sec.*$//')
deps=$(echo "$deps" | sed 's/use m_mrgrnk.*$//')
deps=$(echo "$deps" | sed 's/use m_lookup_table.*$//')
deps=$(echo "$deps" | sed 's/use m_random.*$//')
deps=$(echo "$deps" | sed 's/use m_units_constants.*$//')

# Remove lines without dependencies
deps=$(echo "$deps" | sed 's/^.*: *$//')

# Remove directories
deps=$(echo "$deps" | sed 's/^.*[/]//')

# Fix spacing around ':'
deps=$(echo "$deps" | sed 's/ *: */:/')

# Replace extension
deps=$(echo "$deps" | sed 's/[.]f90/.o/')

# Replace 'use m_xxx' by ' m_xxx.mod'
deps=$(echo "$deps" | sed 's/use \(.*\)$/ \1.mod/')

# Sort lines and remove duplicates
deps=$(echo "$deps" | sort -u)

# Print results
echo "$deps"
