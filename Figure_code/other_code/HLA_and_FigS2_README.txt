## 101723. NB. Gdata package updated yeterday and apparently dropped read.xls function. Thus HLA and other Excel files won't laod with current code. In updated Fig S2, I corrected this in 2 steps, which should be repeated for any script that utilizes read.xls:

1) Use library(readxl) instetad of library(gdata) or require(gdata)
2) use function read_excel
3) Remove stringsAsFactors flag from any excel files previously loaded using read.xls

Check to see that excel files load correctly