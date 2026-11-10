"""Extend a copy of the cohort workbook with all FAM120A follow-up source data."""
from pathlib import Path
import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils.dataframe import dataframe_to_rows

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
OUT=HERE/'publication';FOLLOW=ROOT/'fam120a_followup_20261109'

def main():
    workbook=load_workbook(ROOT/'revision_20261030/publication/Source_Data.xlsx')
    workbook['Readme'].append(['Figure 6 revision','2026-11-10'])
    workbook['Readme'].append(['New panels I-K','FAM120A regions, paired shared targets and LMNB2 original signal'])
    workbook['Readme'].append(['Shared target definition','Both cells mean baseline TPM >= 1; higher/lower TPM <= 2; at least one FAM120A peak per cell'])
    workbook['Readme'].append(['Gene differences','HepG2 minus K562; genes weighted equally; 10000 gene-bootstrap resamples'])
    workbook['Readme'].append(['Signal definition','Absolute value of ENCODE minus-strand signal; no additional normalization; common 0-45 display range in main panel K'])
    workbook['Readme'].append(['Supplements','Supplementary Figure 16: FAM120A sensitivity and locus detail; Supplementary Figure 17: former main-figure boundary/motif examples'])
    source=load_workbook(FOLLOW/'publication/source_data_fam120a.xlsx',read_only=True,data_only=True)
    for sheet in source:
        name='FAM_'+sheet.title
        if name=='FAM_Expression':name='FAM_Baseline_expression'
        assert len(name)<=31 and name not in workbook.sheetnames
        out=workbook.create_sheet(name)
        for row in sheet.iter_rows(values_only=True):out.append(row)
        out.freeze_panes='A2';out.auto_filter.ref=out.dimensions
        for col in out.columns:
            out.column_dimensions[col[0].column_letter].width=min(35,max(12,len(str(col[0].value))+2))
    source.close()
    expression_manifest=pd.read_json(FOLLOW/'source/baseline_expression_manifest.json')
    sheet=workbook.create_sheet('FAM_Expression_sources')
    for row in dataframe_to_rows(expression_manifest,index=False,header=True):sheet.append(row)
    sheet.freeze_panes='A2';sheet.auto_filter.ref=sheet.dimensions
    workbook.save(OUT/'Source_Data_Figure6.xlsx')
    print('Saved',len(workbook.sheetnames),'sheets')

if __name__=='__main__':main()
