# Treatment Gap Radar — code

Full project documentation, methods, findings and run instructions are in the
**[repository README](../README.md)**. Manuscript: **[MANUSCRIPT.md](MANUSCRIPT.md)** ·
slides: **[SLIDES.md](SLIDES.md)** · schema notes: **[PROFILING_NOTES.md](PROFILING_NOTES.md)**.

## Quick start (from this folder)
```bash
conda create -n tgr python=3.12 -y && conda activate tgr
pip install -r ../requirements.txt scikit-learn statsmodels   # + analysis extras
python -m src.pipeline          # build everything (harmonize → indices → rigor → ML → causal)
streamlit run app/dashboard.py  # interactive dashboard
```
Raw Vivli data is **not** in this repo (Data Use Agreement); place it in the parent challenge
folder — paths are in `src/paths.py`. Only aggregate results (n ≥ 30) are published.
