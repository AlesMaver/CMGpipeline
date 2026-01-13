#!/usr/bin/env python3
import requests
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

try:
    response = requests.post("https://mitomap.org/mitomaster/websrvc.cgi", files={"file": open("mtDNA.fasta"),'fileType': ('', 'sequences'),'output':('', 'detail')}, verify=False)
    print(str(response.content, 'utf-8'))
except requests.exceptions.HTTPError as err:
    print("HTTP error: " + err)
