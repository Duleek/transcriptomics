"""
Download files from National Cancer Institute GDC database
"""
import requests
import json
import re

#Use the UUID for the file_id
file_id = "5850111a-db51-4bf0-96ba-fda0423d654a"

data_endpt = "https://api.gdc.cancer.gov/data/{}".format(file_id)

response = requests.get(data_endpt, headers = {"Content-Type": "application/json"})

# The file name can be found in the header within the Content-Disposition key.
response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

with open(file_name, "wb") as output_file:
    output_file.write(response.content)