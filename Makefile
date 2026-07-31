# marineSim — analysis pipeline
#
# The core factorial benchmark takes hours; everything else runs in minutes.
# Case-study targets need the DBCA data (see README, Data availability).

PY      := python3
PYDIR   := python
DATA    ?= data/raw/Bayesiandataset_2025_final.csv

.PHONY: all benchmark ablation verify case figures clean help

help:
	@grep -E '^[a-z]+:.*?##' $(MAKEFILE_LIST) | sed 's/:.*##/\t/'

all: ablation verify case figures  ## everything except the core benchmark

benchmark:  ## core factorial grid (hours)
	cd $(PYDIR) && $(PY) benchmark.py
	cd $(PYDIR) && $(PY) aux_runs.py
	cd $(PYDIR) && $(PY) aggregate.py

ablation:  ## feature-block coalitions and the Shapley decomposition
	cd $(PYDIR) && $(PY) ablation.py 3 0 3 c1
	cd $(PYDIR) && $(PY) ablation.py 3 3 6 c2
	cd $(PYDIR) && $(PY) ablation.py 3 6 9 c3
	cd $(PYDIR) && $(PY) ablation.py 3 9 12 c4
	cd $(PYDIR) && $(PY) ablation.py 3 12 15 c5
	cd $(PYDIR) && $(PY) ablation.py 3 15 18 c6
	cd $(PYDIR) && $(PY) -c "import pandas as pd, glob; \
	  pd.concat([pd.read_csv(f) for f in sorted(glob.glob('../results/ablation_rows_c*.csv'))], \
	  ignore_index=True).to_csv('../results/ablation_rows.csv', index=False)"
	cd $(PYDIR) && $(PY) decompose.py

verify:  ## leakage invariance, GP scaling, prediction-band coverage
	cd $(PYDIR) && $(PY) verify.py

case: $(DATA)  ## case study: LORO, CV regimes, rolling origin, intervals
	cd $(PYDIR) && $(PY) case_study_corrected.py
	cd $(PYDIR) && $(PY) case_study_corrected2.py

$(DATA):
	@echo "Missing $(DATA)."
	@echo "Case-study data are not redistributed; see README, Data availability."
	@exit 1

figures:  ## regenerate figures from stored results
	cd $(PYDIR) && $(PY) render_figs.py

clean:  ## remove intermediate ablation chunks
	rm -f results/ablation_rows_c*.csv
