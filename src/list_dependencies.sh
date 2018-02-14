#!/usr/bin/env bash

# Get all occurrences of use m_... in source files
deps="$(grep -e "use m_" *.f90)"

# Remove Xd files
deps=$(echo "$deps" | sed 's/^\S*_Xd.*$//')

# Remove modules from libraries
deps=$(echo "$deps" | sed 's/use m_a\$D_.*$//')
deps=$(echo "$deps" | sed 's/use m_a2_.*$//')
deps=$(echo "$deps" | sed 's/use m_a3_.*$//')
deps=$(echo "$deps" | sed 's/use m_afivo_.*$//')
deps=$(echo "$deps" | sed 's/use m_particle_core$//')
deps=$(echo "$deps" | sed 's/use m_gas$//')
deps=$(echo "$deps" | sed 's/use m_cross_sec$//')
deps=$(echo "$deps" | sed 's/use m_random$//')
deps=$(echo "$deps" | sed 's/use m_lookup_table$//')
deps=$(echo "$deps" | sed 's/use m_units_constants$//')
deps=$(echo "$deps" | sed 's/use m_.*\$Dd$//')

# Remove comments
deps=$(echo "$deps" | sed 's/!.*$//')

# Remove 'only: ...'
deps=$(echo "$deps" | sed 's/, *only.*$//')

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
