#!/usr/bin/env python3
"""Access DR Comments from HIVdb and arrange them in tsv files.
When masterComments files are updated, mastercomments_version in src/minvar/common.py must be updated too.
"""
import pandas as pd
import io
import re
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Setup the headless browser for an Ubuntu Server environment
options = Options()
options.add_argument('--headless=new') 
options.add_argument('--no-sandbox') 
options.add_argument('--disable-dev-shm-usage') 
options.add_argument('--disable-gpu') 
options.add_argument('--remote-debugging-port=9222')

# Start the invisible browser
driver = webdriver.Chrome(options=options)

classes = ['PI', 'NRTI', 'NNRTI', 'INSTI', 'CAI']  # they are split over NRTI/NNRTI
td = {}
for c in classes:
    print(f"Fetching {c}...")
    driver.get(f'https://hivdb.stanford.edu/dr-summary/comments/{c}/')

    # Wait up to 20 seconds to draw the table
    try:
        WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.TAG_NAME, "table")))
    except:
        print(f"Error: Table did not load for {c} in time.")
        continue

    table = pd.read_html(io.StringIO(driver.page_source))[0]
    table.columns.values[1] = 'Category'
    def parse_condition(val):
        val = str(val) if pd.notnull(val) else ""
        match = re.search(r'(\d*)(\w*)', val)
        return match.group(1) if match else "", match.group(2) if match else ""
    parsed = table['Condition'].apply(parse_condition)
    # extract position/amminoacid
    table['POSITION'], table['AA'] = zip(*parsed)
    td[c] = table

# Shut down the hidden browser
driver.quit()

# merge RTI
td['RTI'] = pd.concat([td['NRTI'], td['NNRTI']])
del td['NRTI']
del td['NNRTI']
# rename INSTI to INI
td['INI'] = td.pop('INSTI')
# rename CAI to CA
td['CA'] = td.pop('CAI')

for k, v in td.items():
    v.to_csv('masterComments_%s.txt' % k, index=False, sep='\t', columns=['POSITION', 'AA', 'Category', 'Comment'])
