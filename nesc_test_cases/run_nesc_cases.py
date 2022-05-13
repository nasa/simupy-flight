import os
import sys
import subprocess

from nesc_testcase_helper import nesc_options

dirname = os.path.dirname(__file__)
if dirname:
    os.chdir(dirname)

savepath = nesc_options["save_relative_path"]
if not os.path.exists(savepath):
    os.makedirs(savepath)

cases_to_run = []
cases_to_run += ["%02d" % num for num in range(1, 12)]
cases_to_run += ["13p%01d" % num for num in range(1, 5)]
cases_to_run += ["%02d" % num for num in range(15, 17)]


def run_case(case_id):
    result = "skipped"
    cmd = ["python", f"nesc_case{case_id}.py"] + sys.argv[1:]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    for line in iter(proc.stdout.readline, ""):
        if "PASSED" in line:
            result = "passed"
        elif "FAILED" in line:
            result = "failed"
        sys.stdout.write(line)
    return result


results = {}
for case_id in cases_to_run:
    print("case number:", case_id, flush=True)
    results[case_id] = run_case(case_id)
    print("\n\n", flush=True)

print("Test summary")
print(22 * "-")
for case_id, result in results.items():
    print(f"case {case_id:8s} {result}")

sys.exit(1 if "failed" in results.values() else 0)
