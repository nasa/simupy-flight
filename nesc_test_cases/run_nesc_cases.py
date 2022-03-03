import os
import sys
from subprocess import call
from nesc_testcase_helper import nesc_options

dirname = os.path.dirname(__file__)
if dirname:
    os.chdir(dirname)

savepath = nesc_options["save_relative_path"]
if not os.path.exists(savepath):
    os.makedirs(savepath)

print("\n\n")
cases_to_run = []
cases_to_run += ["%02d" % num for num in range(1, 12)]
cases_to_run += ["13p%01d" % num for num in range(1, 5)]
cases_to_run += ["%02d" % num for num in range(15, 17)]

for case_id in cases_to_run:
    print("case number:", case_id)
    call(["python", "nesc_case%s.py" % case_id] + sys.argv[1:])

    print("\n\n")
