# marineSim — reproduce the article from the repository root.
# The ablation grid takes hours; everything else runs in minutes to an hour.
RS   := Rscript
A    := analysis
F    := figures/scripts
DATA ?= data/raw/Bayesiandataset_2025_final.csv

.PHONY: help install check verify benchmark decompose gam core si s6 case figures all

help:
	@grep -E '^[a-z]+:.*?##' $(MAKEFILE_LIST) | sed 's/:.*##/\t/'

install:  ## build, check and install the R package
	R CMD build marineSim && _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual marineSim_*.tar.gz && R CMD INSTALL marineSim_*.tar.gz

verify:  ## engine mechanism checks (run first)
	cd $(A) && $(RS) 00_verify_engine.R

benchmark:  ## expanded ablation grid, one core, resumable (~4-5 h)
	cd $(A) && $(RS) 01_ablation_expanded.R

decompose:  ## Shapley tables T1-T5
	cd $(A) && $(RS) 02_decompose.R

gam:  ## spatial GAM on Stage A datasets (~1 h)
	cd $(A) && $(RS) 03_gam_core.R

core:  ## five method core benchmark table (Figure 3, Table 2)
	cd $(A) && $(RS) 04_core_benchmark_table.R

si:  ## S2-S5 data and S6 pooled stRF
	cd $(A) && $(RS) 05_si_figures_data.R && $(RS) 07_s6_pooled.R

s6:  ## sdmTMB comparator (needs sdmTMB; ~2 h)
	cd $(A) && $(RS) 06_s6_sdmtmb.R

case: $(DATA)  ## case study (needs the monitoring data)
	cd $(A) && $(RS) 08_case_study.R

$(DATA):
	@echo "Missing $(DATA). Case study data are not redistributed; see README, Data availability."; exit 1

figures:  ## draw all figures from results/ (ggplot2)
	cd $(F) && $(RS) make_figures_2_4_5.R && $(RS) make_figure_3.R && $(RS) make_figures_S2_S3_S4.R \
	  && $(RS) make_figure_S8.R && $(RS) make_figure_S1.R && $(RS) make_figures_case.R

all: verify benchmark decompose gam core si figures  ## everything except sdmTMB and the case study
