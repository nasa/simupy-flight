import os
from subprocess import call

dirname = os.path.dirname(__file__)
if dirname:
    os.chdir(dirname)

print("\n\n")
for case_num in range(1,11):
    print("case number:", case_num)
    call(["python", "nesc_case%02d.py" % (case_num)])

    print("\n\n")
