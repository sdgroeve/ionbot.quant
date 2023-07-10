# ionbot.quant

Peptide and Protein quantification from ionbot result files with command line [FlashLFQ](https://github.com/smith-chem-wisc/FlashLFQ).

### INSTALL

Install FlashLFQ as a docker image:

    docker pull smithchemwisc/flashlfq:1.0.3

### USAGE

    usage: ionbot2FlashLFQ.py [-h] [-x] folder [folder ...]

    ionbot.quant: conver ionbot result files to FlashLFQ input file

    positional arguments:
        folder      folder with ionbot results

    optional arguments:
        -h, --help  show this help message and exit
        -x          don't filter peptide matches with unexpected modification

Each sample in an experiment is searched by ionbot with all fractions of the same sample in one search (to allow for correct protein inference).

Let's say we have ionbot results for an experiment with 3 samples in folders /home/sample1, /home/sample2 and /home/sample3. 

Then the following script parses the ionbot.first.csv  and ionbot.first.proteins.csv result files for each sample in an experiment:

    python ionbot2FlashLFQ.py /home/sample1 /home/sample2 /home/sample3

It filters for identified PSMs (q-value <= 0.01) not shared between protein groups (read from the ionbot.first.proteins.csv files) and writes a `flashlfq.tsv` result file suitable for quantification by FlashLFQ. 

Peptide matches with an unexpected modification are not used for protein quantification. This can be overruled by setting the `-x` flag.

The `File Name` column in `flashlfq.tsv` points to the spectrum file (can be a fraction) in which the PSM was identified (FlashLFQ in linux reads .mzML files for matching M1 data).

Copy the flashlfq.tsv file to the folder that contains all corresponding .mzML files and run (in that folder):

    docker run --rm -v `pwd`:/mnt/data smithchemwisc/flashlfq:1.0.3 --idt /mnt/data/flashlfq.tsv --rep /mnt/data --mbr

This will quantify al identified peptides and proteins with the Matching Between Runs option enabled.
