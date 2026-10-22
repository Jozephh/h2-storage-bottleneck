# Paste as the first cell of every notebook
import subprocess, hashlib, datetime, pathlib
root = pathlib.Path.cwd().parent
commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=root).stdout.strip()
dirty = subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True, cwd=root).stdout.strip() != ""
inp = hashlib.sha256((root / "data" / "inputs.csv").read_bytes()).hexdigest()[:10]
print(f"Run {datetime.datetime.now():%Y-%m-%d %H:%M} | commit {commit}{' (UNCOMMITTED CHANGES)' if dirty else ''} | inputs.csv {inp}")
