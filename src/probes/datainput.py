import django
import csv
with open('src/cellline_short.csv', newline='') as f:
    csv_reader = csv.DictReader(f)
    datas = [
    (row['CellLineName'], row['PrimarySite'], row['PrimaryHist'])
    for row in csv_reader
   ]
from models import CellLine
for m in datas:
     CellLine(name=m[0], primary_site=m[1], primary_hist=m[2]).save()
