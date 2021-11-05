snakemake -j1 create_participants_files

./rapids -j4 --delete-all-output

./rapids -j26