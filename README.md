# diannReportConvertor
Convert report from DIANN into .msstats and .mzTab
## 1. convert report.tsv into .msstats
`python convert_report2msstats.py convert_report2msstats -ms [mssats] -r [diann_report] -e [exp_design] -u [unimod_csv]`
## 2. convert report.tsv into .mzTab
`python convert_report2mztab.py convert_report2mztab -mz [mztab] -o [openms] -r [diann_report] -pr [report.pr_matrix] -pg [report.pg_matrix] -u [unimod_csv] -un [report.unique_genes_matrix] -f [fasta] -e [experimental_design]`
