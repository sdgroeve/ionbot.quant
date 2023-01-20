# ionbot.quant

Peptide and Protein quantification from ionbot result files with command line [FlashLFQ](https://github.com/smith-chem-wisc/FlashLFQ).

### INSTALL

Install FlashLFQ as a docker image:

    docker pull smithchemwisc/flashlfq:1.0.3

### USAGE

'''
usage: ionbot2FlashLFQ.py [-h] [-x] folder [folder ...]

ionbot.quant: conver ionbot result files to FlashLFQ input file

positional arguments:
  folder      folder with ionbot results

optional arguments:
  -h, --help  show this help message and exit
  -x          don't filter matches with unexpected modification
'''

Each sample in an experiment is searched by ionbot with all fractions of the same sample in one search (to allow for correct protein inference).

Let's say we have ionbot results for an experiment with 3 samples in folders /sample1, /sample2 and /sample3. 
Then the following script parses the ionbot.first.csv (ionbot.first.proteins.csv is required as well to obtain the protein groups) result files for each sample in an experiment:

    python ionbot2FlashLFQ.py sample*/ionbot.first.csv

This script filters for identified PSMs (q-value <= 0.01) not shared between protein groups (read from the ionbot.first.proteins.csv files) and writes a `flashlfq.tsv` result file suitable for quantification by FlashLFQ. The `Peptide Monoisotopic Mass` column is set to the calculated mass for all peptidoforms except those with an unlocalized modification (present in MS1 but not in MS2), where it uses the median precursor mass for each such peptidform over all samples (this is required by FlashLFQ as all `Peptide Monoisotopic Mass` values should be the same for each peptidoform.

The `File Name` column in `flashlfq.tsv` points to the spectrum file (can be a fraction) in which the PSM was identified (FlashLFQ in linux reads .mzML files for matching M1 data).

Copy the flashlfq.tsv file to the folder that contains all corresponding .mzML files and run (in that folder):

    docker run --rm -v `pwd`:/mnt/data smithchemwisc/flashlfq:1.0.3 --idt /mnt/data/flashlfq.tsv --rep /mnt/data --mbr

This will quantify al identified peptides and proteins with the Matching Between Runs option enabled.


### TODO

 - Still some "no quantification" results. 
 - Currently, ionbot.quant only uses the first matches for quantification. Should be easy to add the
   co-eluting.
