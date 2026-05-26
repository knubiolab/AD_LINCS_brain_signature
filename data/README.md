# Data directory

Place large public LINCS/CMAP input files under `data/raw/`. These files are not tracked by git.

Expected layout:

```text
data/raw/
├── GSE92742/
│   ├── GSE92742_Broad_LINCS_sig_info.txt.gz
│   ├── GSE92742_Broad_LINCS_cell_info.txt.gz
│   └── annotated_GSE92742_Broad_LINCS_Level5_COMPZ_n473647x12328.gctx
├── GSE70138/
│   ├── GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz
│   └── GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx
└── CMAP_2020/
    ├── siginfo_beta.txt
    ├── level5_beta_ctl_n58022x12328.gctx
    └── level5_beta_trt_cp_n720216x12328.gctx
```

The scripts accept a custom data directory with `--lincs-dir`.
