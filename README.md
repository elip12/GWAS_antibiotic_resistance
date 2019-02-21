# GWAS_antibiotic_resistance
Methods for conducting genome-wide association studies on microbes to determine genomic origins of antibiotic resistance

## Setting up
1. Clone this repo; go to a folder such as Desktop, and type `git clone https://github.com/elip12/GWAS_antibiotic_resistance.git`
2. Go into the folder. `cd GWAS_antitiotic_resistance`.
3. Make a python virtual environment. I assume you have python3. `python3 -m venv <venv_name>`, where
`<venv_name>` should be replaced with the name you want to give your virtual environment. Mine
is called `v_gwas`.
4. Activate your virtual env. You need to do this every time you quit terminal and restart.
`. <venv_name>/bin/activate`. Note the `.` followed by a space at the beginning.


## Running files
If you want to create the simulated gene pool, change `main.py` and comment out
```
    analyzer = Analyzer(THREADS, K, NUM_SAMPLES, PICKLE_FILE)
    raw = analyzer.load_pool()
    #seqs = analyzer.find_seqs(raw['variants'], SEQS_FILE)
    seqs = analyzer.load(SEQS_FILE)
    analyzer.evaluate_accuracy(seqs, raw)
```
this part, and uncomment out
```
instantiate_pool(PICKLE_FILE)
```
this part. It will take 30 to 90 minutes depending on how shitty your computer is. If you only have
2 cores, I reccommend changing THREADS to 2 (even if you have Intel hyperthreading or equivalent).

If you have a pickle file and want to run the GWAS, switch the commented out stuff above. That
should take 10 to 30 minutes.

You can also tweak anything. Making the simulated data ~30000 bp long will make the simulation and GWAS take close to instant speed.


## Developing

Please do not push to master. Make a new branch and push to it.

