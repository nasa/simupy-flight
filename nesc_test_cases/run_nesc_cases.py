import os
from subprocess import call
from nesc_testcase_helper import parser, save_relative_path

dirname = os.path.dirname(__file__)
if dirname:
    os.chdir(dirname)

if not os.path.exists(save_relative_path):
    os.makedirs(save_relative_path)

args = parser.parse_args()
extra_flags = []
if args.interactive:
    extra_flags.append('--interactive')
if args.simupy_scale:
    extra_flags.append('--simupy-scale')
if args.baseline05:
    extra_flags.append('--baseline05')

print("\n\n")
cases_to_run = ["%02d" % num for num in range(1,12)]
cases_to_run += ["13p%01d" % num for num in range(1,5)]
cases_to_run += ["%02d" % num for num in range(15,17)]

for case_id in cases_to_run:
    print("case number:", case_id)
    call(["python", "nesc_case%s.py" % case_id] + extra_flags)

    print("\n\n")
