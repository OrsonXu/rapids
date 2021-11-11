snakemake -j1 create_participants_files

./rapids -j4 --delete-all-output

./rapids -j26

snakemake -R --dag | dot -Tsvg > dag.svg