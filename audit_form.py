import re

def audit_html_to_pl():
    with open("templates/index.html") as f:
        html = f.read()
    names = re.findall(r'\sname="([^"]+)"', html)
    # filter out standard web ones like sequence_file etc
    names = set(names) - {"sequence_file", "language", "fixed_primer_sequence[]", "fixed_primer_position[]"}

    with open("lava_loop_primer.pl") as f:
        pl = f.read()

    getopts_block = re.search(r"my\s+%optionMap.*?\;", pl, re.DOTALL).group(0)
    options = re.findall(r"\"([a-zA-Z0-9_]+)[=|\"]", getopts_block)
    options = set(options)

    missing = names - options
    if missing:
        print("Found names in HTML NOT in GetOptions:", missing)
    else:
        print("All HTML names are valid Perl options!")

audit_html_to_pl()
