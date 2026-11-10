import re, sys, pathlib
OK = {"HOME","USER","PATH","PWD","SHELL","TMPDIR","HOSTNAME","LANG","RANDOM",
      "BASH_SOURCE","FUNCNAME","PIPESTATUS","IFS","OLDPWD","SHLVL","_",
      "CONDA_PREFIX","CONDA_DEFAULT_ENV","LD_LIBRARY_PATH","CUDA_VISIBLE_DEVICES",
      "REPLY","LINENO","BASHPID","PPID","UID","EUID"}
def strip_comments(src):
    out=[]
    for line in src.splitlines():
        s=line.lstrip()
        if s.startswith("#"):          # whole-line comment (incl. #SBATCH)
            continue
        out.append(line)
    return "\n".join(out)
bad={}
files = sorted(pathlib.Path("slurm").glob("*")) + [pathlib.Path("resume.sh")] \
      + sorted(pathlib.Path("scripts").glob("*.sh"))
for f in files:
    if not f.is_file(): continue
    raw=f.read_text()
    if "set -u" not in raw: continue
    src=strip_comments(raw)
    assigned=set(re.findall(r"^\s*(?:export\s+|local\s+|declare\s+)?([A-Za-z_]\w*)=", src, re.M))
    for m in re.finditer(r"^\s*local\s+([^=\n]+)$", src, re.M):
        assigned |= set(re.findall(r"[A-Za-z_]\w*", m.group(1)))
    for m in re.finditer(r"^\s*local\s+(.+)$", src, re.M):
        assigned |= set(re.findall(r"([A-Za-z_]\w*)=", m.group(1)))
        assigned |= set(re.findall(r"\b([A-Za-z_]\w*)\b(?!=)", m.group(1)))
    assigned |= set(re.findall(r"\bfor\s+([A-Za-z_]\w*)\s+in\b", src))
    assigned |= set(re.findall(r"\bread\s+(?:-r\s+)?([A-Za-z_]\w*)", src))
    # names assigned by an eval of generated KEY="value" lines
    assigned |= set(re.findall(r"print\(f?['\"]\{?(\w+)", src))
    if "eval \"$(" in src:
        assigned |= set(re.findall(r"\b([A-Z][A-Z0-9_]{2,})=", src))
        assigned |= {"ALIGN","PROTEINS","RESULTS","GENOME","RAW"}
    used=set()
    for m in re.finditer(r"\$\{([A-Za-z_]\w*)([^}]*)\}|\$([A-Za-z_]\w*)", src):
        name=m.group(1) or m.group(3); rest=m.group(2) or ""
        if any(t in rest for t in (":-",":=",":?","+","#","%","/")): continue
        used.add(name)
    miss={u for u in used if u not in assigned and u not in OK and not u.startswith("SLURM_")}
    print(f"  {str(f):<34}{'OK' if not miss else 'RISK: '+str(sorted(miss))}")
    if miss: bad[str(f)]=miss
print(f"\nfiles with risk: {len(bad)}")
sys.exit(1 if bad else 0)
