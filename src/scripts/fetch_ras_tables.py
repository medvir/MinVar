#!/usr/bin/env python3
"""Parse tables from http://onlinelibrary.wiley.com/doi/10.1002/hep.27934/full and write them to files."""
import shlex
import os.path
import subprocess
from bs4 import BeautifulSoup

# Remove page.html to download page again

s = 'wget -O page.html http://onlinelibrary.wiley.com/doi/10.1002/hep.27934/full'

if not os.path.exists('page.html'):
    subprocess.call(shlex.split(s))


def get_table(number):
    """Workhorse."""
    with open("page.html") as fp:
        soup = BeautifulSoup(fp)

    tables = soup.find_all("table",
                           attrs={"class": "table table--article-section"})
    table = tables[number - 1]
    head = [th.get_text() for th in table.find_all("tr")[1].find_all("th")]
    datasets = []
    for row in table.find_all("tr")[3:]:
        ds = [td.get_text().replace(u'\xa0', '') for td in row.find_all("td")]
        datasets.append(ds)
    with open("RAS_table%d.csv" % number, 'w') as oh:
        oh.write(','.join(head) + '\n')
        for r in datasets:
            oh.write(','.join(r) + '\n')


get_table(1)
get_table(2)
get_table(3)
