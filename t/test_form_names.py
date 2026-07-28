import re
import sys

def run_test():
    with open("templates/index.html") as f:
        html = f.read()

    # Extract all name="..."
    names = set(re.findall(r'\sname="([^"]+)"', html))
    
    # These are handled separately by flask (not passed directly as simple kwargs to perl)
    ui_fields = {
        "sequence_file", "language", "csrf_token", "output_name", "fasta_file", 
        "lamp_mode", "script_type", "fixed_primer_strict",
        "fp_type[]", "fp_seq[]", "fp_pos[]", "fixed_primer_sequence[]", "fixed_primer_position[]"
    }
    
    names = names - ui_fields

    # Extract valid parameters from lava_flask_app.py
    with open("lava_flask_app.py") as f:
        flask_app = f.read()
        
    all_valid_params = set()
    for var_name in ["common_params", "loop_only_params", "stem_only_params"]:
        match = re.search(var_name + r"\s*=\s*\{(.*?)\}", flask_app, re.DOTALL)
        if match:
            valid_params_str = match.group(1)
            # Extract string literals
            valid_params = set(re.findall(r"['\"]([^'\"]+)['\"]", valid_params_str))
            all_valid_params.update(valid_params)
    
    # Also grab keys from param_mapping if missing
    mapping_match = re.search(r"param_mapping\s*=\s*\{(.*?)\}", flask_app, re.DOTALL)
    if mapping_match:
        keys = re.findall(r"['\"]([^'\"]+)['\"]\s*:\s*['\"]", mapping_match.group(1))
        all_valid_params.update(keys)

    # Check if there are any names in HTML not in valid_params
    missing = names - all_valid_params
    
    if missing:
        print("FAIL: Found fields in index.html that are NOT in valid Perl param lists (common_params, etc.):")
        for m in missing:
            print(f"  - {m}")
        sys.exit(1)
    else:
        print("PASS: All HTML fields are valid Perl params.")
        sys.exit(0)

if __name__ == "__main__":
    run_test()
