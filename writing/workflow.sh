#!/bin/bash

# For use in compiling the after_peer_review files and tracking the changes
# from the before_peer_review files

# Note - ```latexmk``` detects targets and compiles everything
# Note - ```latexdiff``` tracks changes
# https://www.overleaf.com/learn/latex/Articles/Using_Latexdiff_For_Marking_Changes_To_Tex_Documents

prev_source=before_peer_review
new_source=rev_2
diff_folder=./diff_${new_source}_cfont_mod

cd ./${new_source}
latexmk
rm *.fdb_latexmk # work-around hack
latexmk
cd ../

rm -r $diff_folder
mkdir -p ./${diff_folder}
rsync -rav --delete ./${new_source}/images/ ./${diff_folder}/images/

for file in "elsarticle.cls" "elsarticle-num.bst"; do
    rsync ./${new_source}/${file} ./${diff_folder}/
done

for file in "main.bib" "main.bbl"; do
latexdiff ./${prev_source}/${file} ./${new_source}/${file} > ./${diff_folder}/${file}
done

for file in main supp figures highlights; do
latexdiff --flatten --packages=xr,hyperref,biblatex ./${prev_source}/${file}.tex ./${new_source}/${file}.tex > ./${diff_folder}/${file}.tex
done

# latexdiff options:  --type=CTRADITIONAL, for example
# manually change the title entries before running latexmk



