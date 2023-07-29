production (previously rev_2) contains the most up to date revision

diff_rev_2 folder contains the output of running latexdiff
(comparing rev_2 with before_peer_review) as shown in the workflow script

In highlighting
    After running workflow (for latexdiff), I manually removed latexdiff formatting from the titles
    (which was appearing to cause problems)

    I substituted
    \providecommand{\DIFdel}[1]{} %Don't show deleted text
    in place for whatever the previous
    \providecommand{\DIFdel}<something else>
    preamble was

    I ran latexmk (note this ran to completion without errors, but cross references didn't work)

    I manually hard-coded cross reference numbers in the main .tex file

    I re-ran latexmk

    Cross references did work in supp and figure files though

v2 was built from source (workflow_source), as it should be, and as it works locally
v3 - manually hard coded references to all figures and supplementary material
    (brute force solution)

FOR MY FUTURE REFERENCE_________________________________________________________

latexmk is apparently imperfect - when I tried compiling my latex document from
only source files, cross references didn't work. When I tried compiling from
source files + .aux files, it did work. The strange thing though is that the .aux
files that latexmk builds from scratch are the same as those that are needed for
successful compilation. So as a work-around, isolate the source files, run
latexmk, delete the *.fdb_latexmk files, and re-run latexmk

Also, Elsevier submissions don't allow latex macros (\newcommand)
