# afivo-pic
Parallel PIC simulations of electric discharges, with adaptive mesh refinement

# Getting started

The first time, run these commands to download `afivo`, `particle_core` etc.:

    git submodule init
    git submodule update

Then compile the examples in `programs` by typing `make` in the example folders. Afterwards, the examples can be run with:

    ./apic my_config_file.cfg

# Updating the code

First get the latest version:

    git pull
    git submodule update

And then type `make' in the respective program to compile again.

