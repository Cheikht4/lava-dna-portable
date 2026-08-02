import requests
import re

url = "http://127.0.0.1:5000"

s = requests.Session()
r = s.get(url)
print("GET status:", r.status_code)

match = re.search(r'name="csrf_token" value="([^"]+)"', r.text)
if match:
    csrf_token = match.group(1)
    print("Found CSRF Token:", csrf_token)
else:
    print("No CSRF token found in index.html!")
    exit(1)

with open("test_upload.fasta", "w") as f:
    f.write(">seq1\nACGT\n>seq2\nACGT\n")

files = {'fasta_file': open("test_upload.fasta", 'rb')}
data = {'csrf_token': csrf_token}

r2 = s.post(url + "/upload", files=files, data=data)
print("POST status:", r2.status_code)
if r2.status_code == 400:
    print("Response text:", r2.text)
