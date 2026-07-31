# Releasing marineSim to GitHub and Zenodo

A checklist for archiving this repository so the article can cite a permanent DOI.

**Why both?** Methods in Ecology and Evolution is a British Ecological Society journal.
BES policy does not accept archiving on GitHub alone, because GitHub does not mint DOIs
and its contents can be changed or deleted. Zenodo mints a permanent DOI and freezes the
contents. GitHub is where the code lives; Zenodo is what the paper cites.

Order matters throughout. Step 2.2 in particular must happen **before** step 2.3, or you
will get no DOI and will have to redo the release.

---

## Step 0 — Before you push anything

Three of these are blockers. Do not skip them.

### 0.1 Decide what you are legally allowed to publish  **[BLOCKER]**

The case study uses long term *Posidonia sinuosa* monitoring data (17,935 quadrats,
five marine parks). That data was collected by monitoring programmes and is **not
automatically yours to redistribute**. Before it goes anywhere public:

- Confirm in writing with the data custodians (the monitoring programme leads and the
  relevant state agency) whether the site-year aggregate may be released, and under what
  licence.
- If yes: add it under `data/` with its own licence file, and remove the exclusion in
  `.gitignore`.
- If no: keep it excluded. The code still archives fine; your Data Availability Statement
  must then say plainly that the case-study data are available from the custodians on
  request, and name the contact point. Do not write "available on request" and leave it
  at that, and do not quietly ship the data hoping nobody checks.

The simulation is fully synthetic and has no such constraint. It regenerates from seed.

### 0.2 Run the R package check  **[BLOCKER]**

The `marineSim` package has never been executed. It compiles as written, but that is not
the same as working.

```r
install.packages(c("devtools", "ranger", "mgcv", "FNN"))
devtools::check("marineSim")
devtools::test("marineSim")
```

Fix whatever it reports. Archiving a package that fails `R CMD check` is worse than
archiving no package, because a reviewer will run it.

### 0.3 Confirm the outputs are actually committed

`.gitignore` previously excluded `figures/*.png` and `results/*.{csv,json}`, which would
have pushed empty directories. That is fixed, but verify after your first commit:

```bash
git ls-files figures/ | wc -l     # expect 20
git ls-files results/ | wc -l     # expect 25
```

If those come back as 0, your `.gitignore` is still excluding them.

### 0.4 Fill in the placeholders

| File | Placeholder | Replace with |
|---|---|---|
| `CITATION.cff` | `<GITHUB-USER>` | your GitHub account name |
| `CITATION.cff` | `date-released` | the date you cut the release |
| `.zenodo.json` | `GITHUB-USER` | your GitHub account name |
| `README.md` | Citation section | the Zenodo DOI, once you have it (step 2.4) |
| Manuscript | `10.xxxx/figshare.placeholder` | the Zenodo DOI (step 3) |
| Manuscript | `github.com/eafrifayamoah/marineSim` | your real repository URL |

Co-author ORCIDs are marked `TODO` in `CITATION.cff`. Adding them means Zenodo credits
every author properly rather than just you.

---

## Step 1 — GitHub

### 1.1 Create the repository

On github.com: **New repository** → name it `marineSim` → **Public** (Zenodo cannot see
private repositories) → do **not** initialise with a README, licence or .gitignore, since
this repository already has them.

### 1.2 Push

From the repository root:

```bash
git init
git branch -M main
git add .
git commit -m "marineSim: stRF simulation benchmark, R package and seagrass case study"
git remote add origin https://github.com/<GITHUB-USER>/marineSim.git
git push -u origin main
```

### 1.3 Check the result in a browser

- The README renders.
- `figures/` and `results/` are populated, not empty.
- GitHub shows a **Cite this repository** button on the right (this means it parsed
  `CITATION.cff`; if it is missing, the YAML has a syntax error).
- No data file you are not permitted to share has been pushed. If one slipped through,
  deleting it in a later commit is **not enough** — it stays in the history. You would
  need to rewrite history or, more safely, delete the repository and start again.

---

## Step 2 — Zenodo

### 2.1 Log in and link ORCID

Go to zenodo.org → **Log in** → choose **Log in with GitHub** → authorise.

In Zenodo, open your profile and connect your ORCID (`0000-0003-1741-9249`). This makes
the deposition appear on your ORCID record automatically.

### 2.2 Switch the repository on  **[DO THIS BEFORE THE RELEASE]**

In Zenodo: your username → **GitHub** → find `marineSim` in the list → flip the toggle to
**On**.

This is the step people get wrong. **Zenodo only archives releases created after the
toggle is on.** If you create the GitHub release first and switch Zenodo on afterwards,
nothing is archived and no DOI is minted; you have to delete the release and cut a new
one. If your repository is not listed, click **Sync now** — Zenodo only lists public
repositories.

### 2.3 Cut the release on GitHub

Back on GitHub: **Releases** → **Create a new release**.

- **Tag**: `v0.3.0` (match `Version:` in `marineSim/DESCRIPTION`, and keep the `v`)
- **Title**: `marineSim v0.3.0 — stRF for marine SDM under data limitation`
- **Description**: one paragraph on what the release contains, and note that it
  accompanies the MEE submission.
- **Publish release.**

Zenodo receives the webhook, ingests the snapshot and mints the DOI. Allow a couple of
minutes.

### 2.4 Collect the DOIs — there are two, and the difference matters

Open the new record on Zenodo. You will see:

- a **concept DOI** ("Cite all versions"), which always resolves to the newest version;
- a **version DOI**, which resolves to this exact snapshot and nothing else.

**Cite the version DOI in the paper.** Reproducibility is the point: a reader must land
on precisely the code that produced Table 3, not on whatever the code became two years
later. Put the concept DOI in the README badge, where "latest" is what you want.

### 2.5 Tidy the Zenodo metadata

`.zenodo.json` pre-fills most of it, but check on the record's **Edit** page:

- Upload type is **Software**.
- Licence is **MIT** (Zenodo defaults to CC-BY, which is a documents licence, not a code
  licence — change it if it did not pick up the file).
- All five authors are listed with affiliations and ORCIDs.
- Add the article DOI under **Related identifiers** as *is supplement to* once the paper
  is accepted.

Metadata is editable forever. **The files are not.** A published Zenodo record cannot have
its files swapped — a correction means a new version, and a new version DOI.

---

## Step 3 — Wire the DOI back into the manuscript

Two placeholders currently sit in the manuscript and must both change:

1. `10.xxxx/figshare.placeholder` → your Zenodo **version** DOI. Note the manuscript text
   currently says *figshare*; that word has to change to *Zenodo* as well, not just the
   number.
2. `github.com/eafrifayamoah/marineSim` → your real repository URL.

Suggested Data Availability Statement, assuming the seagrass data cannot be redistributed
(adjust if step 0.1 goes the other way):

> All code required to reproduce the simulation study, the case study analyses and every
> figure and table is archived on Zenodo (doi:10.5281/zenodo.XXXXXXX) and developed at
> https://github.com/<GITHUB-USER>/marineSim. The simulation study is fully synthetic and
> regenerates from the seeds recorded in the archive. The *Posidonia sinuosa* monitoring
> data underpinning the case study are held by [custodian] and are available on request
> from [contact]; they are not redistributed here because they are third-party monitoring
> data.

Add the badge to the top of `README.md` (concept DOI, so it tracks the latest version):

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

---

## Step 4 — Revisions

Reviewers will ask for changes, and the code will change with them.

1. Commit and push to `main` as usual.
2. Bump `Version:` in `marineSim/DESCRIPTION` (e.g. `0.3.1`) and `version:` in
   `CITATION.cff`.
3. Cut a new GitHub release (`v0.3.1`).
4. Zenodo automatically archives it as a **new version** of the same record.

The concept DOI is unchanged. A **new version DOI** is minted, and that is the one the
final accepted manuscript should cite. Update the DOI in the paper at proof stage so the
archived code matches the published text.

---

## Quick reference

| | Purpose | Which DOI |
|---|---|---|
| GitHub | Development, issues, visibility | none |
| Zenodo | Permanent archive, what MEE requires | version DOI in the paper |
| Zenodo concept DOI | Always points at latest | README badge |

**The one that bites:** switch Zenodo on (2.2) *before* cutting the release (2.3).
