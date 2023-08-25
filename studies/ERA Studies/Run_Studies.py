import subprocess

files = ["ERA_Study_1.py", "ERA_Study_2.py", "ERA_Study_3.py", "ERA_Study_4.py", "ERA_Study_5.py"]

for file in files:
    subprocess.run(["python", file])