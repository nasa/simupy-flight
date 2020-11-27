import os
from subprocess import call
from nesc_testcase_helper import parser

dirname = os.path.dirname(__file__)
if dirname:
    os.chdir(dirname)

args = parser.parse_args()
extra_flags = []
if args.interactive:
    extra_flags.append('--interactive')
if args.simupy_scale:
    extra_flags.append('--simupy_scale')

print("\n\n")
cases_to_run = ["%02d" % num for num in range(1,11)]

for case_id in cases_to_run:
    print("case number:", case_id)
    call(["python", "nesc_case%s.py" % case_id] + extra_flags)

    print("\n\n")
